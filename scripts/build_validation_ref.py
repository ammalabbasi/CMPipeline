#!/usr/bin/env python3
"""Alignment validation, step 2: build the competitive reference (public genomes only).

For every species in validation_species.tsv, take RefSeq genomes by taxid (species_taxid or
taxid; falling back to an exact organism-name match when taxonomy has moved since the Kraken
DB was built):
  * latest assemblies only; Complete Genome if the species has any, else its best available
    level (Chromosome > Scaffold > Contig); RefSeq reference genome first, then a seeded shuffle;
    at most --max_genomes per species
  * downloaded from NCBI with retries and a full gzip integrity check
  * de-duplicated within species with skani: a genome is dropped if it is >= --dedup_ani to one
    already kept (so the reference genome, first in line, is always kept)
  * ANI-curated per genus (validation_ani.py): genomes filed under the wrong species are moved
    to the species they belong to, and species whose anchors meet at >= --group_ani form an ANI
    group, the scorer's fallback level
Human T2T-CHM13v2.0 is added to the same reference so leftover human reads compete.

Outputs (in --out_dir):
  validation_reference.fna     all contigs, named <label>|<accession>|<contig>
  contig_labels.tsv            contig, label, accession, length, ani_group ('' when ungrouped)
  validation_reference_manifest.tsv   one row per genome: label, species, taxid, accession,
                               assembly_level, refseq_category, md5, url, source summary date,
                               filed_label, ani_action, ani_score_own, ani_best_label, ani_score_best
  missing_species.tsv          species that got no genome, and why
The run fails if more than --max_missing_frac of species got no genome.
"""
import argparse, collections, concurrent.futures as cf, csv, datetime as dt, gzip, hashlib, os, random, shutil
import subprocess, sys, tempfile, time, urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validation_common import HUMAN_LABEL
import validation_ani

NCBI = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq"
GROUPS = ("bacteria", "archaea", "fungi", "protozoa")
T2T_URL = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/"
           "GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz")
LEVEL_RANK = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}


def fetch(url, dest, tries=4):
    """Download url -> dest atomically, verifying gzip integrity for .gz files."""
    tmp = dest + ".part"
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as out:
                shutil.copyfileobj(r, out, 1 << 20)
            if dest.endswith(".gz"):
                with gzip.open(tmp, "rb") as fh:
                    while fh.read(1 << 22):
                        pass
            os.replace(tmp, dest)
            return True
        except Exception as e:                       # noqa: BLE001 - network: retry anything
            err = e
            time.sleep(5 * (i + 1))
    if os.path.exists(tmp):
        os.remove(tmp)
    print(f"WARN: download failed after {tries} tries: {url} ({err})", file=sys.stderr)
    return False


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def load_summaries(cache, max_age_days, snapshot=None):
    """RefSeq assembly summaries, cached under cache/; refreshed when older than max_age_days.

    `snapshot` (YYYY-MM) is part of the reference's cache key (decision 26, audit B04). A summary
    from another month is refreshed; if the snapshot names a past month, it cannot be rebuilt
    (NCBI serves only the current summary), so the build stops instead of mislabelling a newer panel.
    """
    os.makedirs(cache, exist_ok=True)
    this_month = dt.date.today().strftime("%Y-%m")
    if snapshot and snapshot != this_month:
        sys.exit(f"ERROR: --snapshot {snapshot} is not the current month ({this_month}), and a past RefSeq "
                 "snapshot cannot be re-downloaded. Its reference is not in the cache; pin a saved "
                 "reference with --validation_ref, or drop --validation_refseq_snapshot to build this month's.")
    by_taxid, by_name, dates = {}, {}, {}
    for g in GROUPS:
        path = os.path.join(cache, f"assembly_summary_{g}.txt")
        stale = not os.path.exists(path) or \
            (time.time() - os.path.getmtime(path)) > max_age_days * 86400 or \
            (snapshot and dt.date.fromtimestamp(os.path.getmtime(path)).strftime("%Y-%m") != snapshot)
        if stale and not fetch(f"{NCBI}/{g}/assembly_summary.txt", path):
            if not os.path.exists(path):
                sys.exit(f"ERROR: could not download the RefSeq {g} assembly summary")
        dates[g] = dt.date.fromtimestamp(os.path.getmtime(path)).isoformat()
        if snapshot and not dates[g].startswith(snapshot):
            sys.exit(f"ERROR: the RefSeq {g} summary is from {dates[g]}, not snapshot {snapshot} "
                     "(download failed?); refusing to build a mislabelled reference")
        with open(path) as fh:
            fh.readline()
            header = fh.readline().lstrip("#").rstrip("\n").split("\t")
            for r in csv.DictReader(fh, fieldnames=header, delimiter="\t"):
                if r["version_status"] != "latest" or r["ftp_path"] in ("", "na"):
                    continue
                r["_group"] = g
                for t in {r["species_taxid"], r["taxid"]}:
                    by_taxid.setdefault(t, []).append(r)
                by_name.setdefault(" ".join(r["organism_name"].split()[:2]), []).append(r)
    return by_taxid, by_name, dates


def choose(rows, max_genomes, seed):
    """Best available level only; reference genome first, then a seeded shuffle."""
    rows = list({r["assembly_accession"]: r for r in rows}.values())
    best = min(LEVEL_RANK.get(r["assembly_level"], 9) for r in rows)
    rows = [r for r in rows if LEVEL_RANK.get(r["assembly_level"], 9) == best]
    ref = [r for r in rows if r["refseq_category"] == "reference genome"]
    rest = sorted((r for r in rows if r["refseq_category"] != "reference genome"),
                  key=lambda r: r["assembly_accession"])
    random.Random(seed).shuffle(rest)
    return (ref + rest)[:max_genomes]


def dedup(label, files, skani, ani, threads, tmp):
    """Greedy de-duplication in order; keeps the first genome of every >= ani cluster."""
    if len(files) < 2:
        return files
    lst = os.path.join(tmp, label + ".list")
    with open(lst, "w") as fh:
        fh.write("\n".join(files) + "\n")
    out = subprocess.run([skani, "triangle", "-l", lst, "-E", "-t", str(threads)],
                         capture_output=True, text=True, check=True).stdout.splitlines()
    pair = {}
    for line in out[1:]:
        p = line.split("\t")
        pair[(p[0], p[1])] = pair[(p[1], p[0])] = float(p[2])
    kept = []
    for f in files:
        if all(pair.get((f, k), 0.0) < ani for k in kept):
            kept.append(f)
    return kept


def skani_pairs(files, skani, threads, tmp, name):
    """skani triangle -E over files -> {(file_a, file_b): ani}, both orders."""
    if len(files) < 2:
        return {}
    lst = os.path.join(tmp, name + ".list")
    with open(lst, "w") as fh:
        fh.write("\n".join(files) + "\n")
    out = subprocess.run([skani, "triangle", "-l", lst, "-E", "-t", str(threads)],
                         capture_output=True, text=True, check=True).stdout.splitlines()
    pair = {}
    for line in out[1:]:
        p = line.split("\t")
        pair[(p[0], p[1])] = pair[(p[1], p[0])] = float(p[2])
    return pair


def ani_curate(kept_rows, skani, threads, tmp, anchor_ani, margin, group_ani):
    """Relabel and group per genus (validation_ani.py). -> (rows with their final label, {label: group},
    {accession: relabel record}). Dropped genomes are left out of the returned rows."""
    by_genus = collections.defaultdict(list)
    for item in kept_rows:
        by_genus[validation_ani.genus_of(item[0]["label"])].append(item)
    final, group_of, record = [], {}, {}
    for genus, items in sorted(by_genus.items()):
        if len({i[0]["label"] for i in items}) < 2:
            final += [(i[0]["label"],) + i for i in items]
            continue
        pair = skani_pairs([i[3] for i in items], skani, threads, tmp, "genus_" + genus)
        acc_of = {i[3]: i[1]["assembly_accession"] for i in items}
        ani = {(acc_of[a], acc_of[b]): v for (a, b), v in pair.items()}
        genomes = [(i[1]["assembly_accession"], i[0]["label"]) for i in items]
        for r in validation_ani.relabel(genomes, ani, anchor_ani, margin):
            record[r["accession"]] = r
        moved = []
        for i in items:
            r = record[i[1]["assembly_accession"]]
            if r["label"]:
                moved.append((r["label"],) + i)
        final += moved
        taxid_of = {i[0]["label"]: i[0]["taxid"] for i in items}
        group_of.update(validation_ani.groups([(m[2]["assembly_accession"], m[0]) for m in moved],
                                              ani, taxid_of, anchor_ani, group_ani))
    acts = collections.Counter(r["action"] for r in record.values())
    print("ANI curation: " + (", ".join(f"{k}={v}" for k, v in sorted(acts.items())) or "no multi-species genus")
          + f"; {len(set(group_of.values()))} ANI group(s) at {group_ani:g}%")
    return final, group_of, record


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True)
    ap.add_argument("--cache_dir", required=True, help="shared cache for summaries, genomes and T2T")
    ap.add_argument("--out_dir", default=".")
    ap.add_argument("--max_genomes", type=int, default=50)
    ap.add_argument("--dedup_ani", type=float, default=99.5)
    ap.add_argument("--max_missing_frac", type=float, default=0.10)
    ap.add_argument("--summary_max_age_days", type=int, default=30)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--skani", default="skani")
    ap.add_argument("--seed", type=int, default=20260923)
    ap.add_argument("--anchor_ani", type=float, default=95.0,
                    help="ANI that defines a species' anchor cluster and the bar for relabelling")
    ap.add_argument("--relabel_margin", type=float, default=1.0,
                    help="a genome moves only if its new species beats its filed one by this much ANI")
    ap.add_argument("--group_ani", type=float, default=92.0,
                    help="species whose anchors meet at this ANI form one fallback group (0 = none)")
    ap.add_argument("--max_ref_gb", type=float, default=40,
                    help="refuse a planned reference larger than this many Gb of sequence (audit B05)")
    ap.add_argument("--snapshot", default=None, help="RefSeq snapshot month (YYYY-MM) in the cache key")
    a = ap.parse_args(argv)

    species = list(csv.DictReader(open(a.species), delimiter="\t"))
    by_taxid, by_name, dates = load_summaries(a.cache_dir, a.summary_max_age_days, a.snapshot)
    gdir = os.path.join(a.cache_dir, "genomes"); os.makedirs(gdir, exist_ok=True)

    plan, missing = [], []
    for s in species:
        rows = by_taxid.get(s["taxid"], []) if s["taxid"] != "NA" else []
        how = "taxid"
        if not rows:
            rows, how = by_name.get(s["species"], []), "name"
        if not rows:
            missing.append((s, "no RefSeq assembly for this taxid or name")); continue
        for r in choose(rows, a.max_genomes, a.seed):
            leaf = r["ftp_path"].rstrip("/").rsplit("/", 1)[1]
            url = r["ftp_path"].rstrip("/").replace("ftp://", "https://") + f"/{leaf}_genomic.fna.gz"
            plan.append((s, r, url, os.path.join(gdir, leaf + "_genomic.fna.gz"), how))

    # Size the reference BEFORE downloading or indexing anything (audit B05). RefSeq summaries carry
    # genome_size; a genome without it is counted at 5 Mb (a typical bacterium). Plus T2T (3.1 Gb).
    known = [int(r["genome_size"]) for _, r, _, _, _ in plan if (r.get("genome_size") or "").isdigit()]
    planned_bases = sum(known) + 5_000_000 * (len(plan) - len(known)) + 3_117_275_501
    print(f"planned reference: {planned_bases / 1e9:.1f} Gb before de-duplication "
          f"({len(plan) - len(known)} genome(s) without a size, counted at 5 Mb); cap {a.max_ref_gb:g} Gb")
    if planned_bases > a.max_ref_gb * 1e9:
        sys.exit(f"ERROR: the planned validation reference is {planned_bases / 1e9:.1f} Gb, above "
                 f"--validation_max_ref_gb {a.max_ref_gb:g}. Lower --validation_max_genomes, raise "
                 "--validation_min_support (fewer species), or raise the cap (index memory grows with it).")

    todo = [(url, dest) for _, _, url, dest, _ in plan if not os.path.exists(dest)]
    print(f"{len(species)} species, {len(plan)} genomes planned, {len(todo)} to download")
    with cf.ThreadPoolExecutor(max_workers=a.threads) as ex:
        ok = dict(zip([d for _, d in todo], ex.map(lambda t: fetch(*t), todo)))

    t2t = os.path.join(a.cache_dir, "T2T-CHM13v2.0.fna.gz")
    if not os.path.exists(t2t) and not fetch(T2T_URL, t2t):
        sys.exit("ERROR: could not download the human T2T reference")

    tmp = tempfile.mkdtemp(prefix="valref_")
    kept_rows = []
    by_label = {}
    for s, r, url, dest, how in plan:
        if os.path.exists(dest):
            by_label.setdefault(s["label"], []).append((s, r, url, dest, how))
    for s in species:
        if s["label"] not in by_label and not any(m[0] is s for m in missing):
            missing.append((s, "all genome downloads failed"))
    for label, items in by_label.items():
        files = dedup(label, [i[3] for i in items], a.skani, a.dedup_ani, a.threads, tmp)
        kept_rows += [i for i in items if i[3] in files]
    final_rows, group_of, ani_record = ani_curate(kept_rows, a.skani, a.threads, tmp, a.anchor_ani,
                                                  a.relabel_margin, a.group_ani)
    if a.group_ani <= 0:
        group_of = {}
    shutil.rmtree(tmp)
    # An unnamed "sp." species whose every genome moved to a named one leaves the reference.
    still = {r[0] for r in final_rows}
    for s in species:
        if s["label"] in by_label and s["label"] not in still:
            to = sorted({rec["label"] for rec in ani_record.values() if rec["filed_label"] == s["label"]} - {""})
            missing.append((s, "every genome relabelled by ANI to " + (", ".join(to) or "none (dropped)")))

    os.makedirs(a.out_dir, exist_ok=True)
    fna = open(os.path.join(a.out_dir, "validation_reference.fna"), "w")
    labels = open(os.path.join(a.out_dir, "contig_labels.tsv"), "w")
    labels.write("contig\tlabel\taccession\tlength\tani_group\n")

    def emit(path, label, acc):
        name, n = None, 0
        grp = group_of.get(label, "")
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    if name:
                        labels.write(f"{name}\t{label}\t{acc}\t{n}\t{grp}\n")
                    name, n = f"{label}|{acc}|{line[1:].split()[0]}", 0
                    fna.write(f">{name}\n")
                else:
                    n += len(line.strip()); fna.write(line)
        if name:
            labels.write(f"{name}\t{label}\t{acc}\t{n}\t{grp}\n")

    with open(os.path.join(a.out_dir, "validation_reference_manifest.tsv"), "w") as man:
        # label = the label the genome is scored under; filed_label = RefSeq's. Dropped genomes
        # (drop_ambig) are listed with an empty label and are not in the reference.
        man.write("label\tspecies\ttaxid\tmatched_by\taccession\tassembly_level\trefseq_category\t"
                  "md5\turl\tsummary_date\tfiled_label\tani_action\tani_score_own\tani_best_label\t"
                  "ani_score_best\n")
        final_label = {r[2]["assembly_accession"]: r[0] for r in final_rows}
        for s, r, url, dest, how in kept_rows:
            acc = r["assembly_accession"]
            rec = ani_record.get(acc)
            label = final_label.get(acc, "")
            if label:
                emit(dest, label, acc)
            ani_cols = (f"{rec['action']}\t{rec['score_own']:.2f}\t{rec['best_label']}\t{rec['score_best']:.2f}"
                        if rec else "NA\tNA\tNA\tNA")
            man.write(f"{label}\t{s['species']}\t{s['taxid']}\t{how}\t{acc}\t"
                      f"{r['assembly_level']}\t{r['refseq_category']}\t{md5(dest)}\t{url}\t{dates[r['_group']]}\t"
                      f"{s['label']}\t{ani_cols}\n")
        emit(t2t, HUMAN_LABEL, "GCF_009914755.1")
        man.write(f"{HUMAN_LABEL}\tHomo sapiens\t9606\tfixed\tGCF_009914755.1\tComplete Genome\t"
                  f"reference genome\t{md5(t2t)}\t{T2T_URL}\tNA\t{HUMAN_LABEL}\tNA\tNA\tNA\tNA\n")
    fna.close(); labels.close()

    with open(os.path.join(a.out_dir, "missing_species.tsv"), "w") as out:
        out.write("label\tspecies\ttaxid\treason\n")
        for s, why in missing:
            out.write(f"{s['label']}\t{s['species']}\t{s['taxid']}\t{why}\n")
    n_genomes = len(final_rows)
    # labels that left by ANI relabelling are listed but are not missing genomes
    n_missing = sum(1 for _, why in missing if not why.startswith("every genome relabelled"))
    frac = n_missing / len(species) if species else 0
    print(f"reference: {n_genomes} genomes for {len(by_label)} species (+ human T2T); "
          f"{n_missing} species missing ({frac:.1%}); {len(plan) - n_genomes} dropped as duplicates or failed")
    if frac > a.max_missing_frac:
        sys.exit(f"ERROR: {n_missing} of {len(species)} species ({frac:.1%}) got no genome, above "
                 f"--validation_max_missing_frac {a.max_missing_frac:g}; see missing_species.tsv")


if __name__ == "__main__":
    main()
