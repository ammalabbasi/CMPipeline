"""ANI curation of the validation reference (decision 36; ANI check of 6 Oct 2026, decided 26 Oct 2026).

RefSeq species labels are not always right. In the CRC50 reference, 4 genomes filed as
Fusobacterium nucleatum are F. animalis or F. vincentii (97-98% median ANI to those species,
about 91% to nucleatum). Every read of a true F. animalis then ties between two labels, so it is
not unique, and F. animalis cells came out `insufficient` (39% unique reads; 85% after relabelling,
public-genome simulation, job 13670581). PRISM, by comparison, uses the filed NCBI labels as they are.

Two steps, both per genus, on the genomes kept after the 99.5% de-duplication:

relabel(...)   1. Anchor of label L = the >= anchor_ani single-linkage cluster of L's own genomes
                  that contains L's first genome (the RefSeq reference genome when there is one).
               2. score(g, M) = median ANI of genome g to M's anchor genomes (g excluded; pairs skani
                  did not report count as 0).
               3. keep        the best label is L (or no label reaches anchor_ani: `divergent`, or
                              `keep_singleton` when L has one genome; kept under L, flagged)
                  relabel     the best label M != L, score >= anchor_ani, and beats score(g, L) by >= margin
                  drop_ambig  the best label M != L reaches anchor_ani but is within margin of L: dropped
                  Unnamed "_sp._" labels are never targets. A NAMED label's seed never moves; an unnamed
                  label's may, and the label then leaves the reference (its reads were a named
                  species': 4 Bacteroides and 4 Escherichia "sp." genomes in CRC50).
               The unnamed 18-genome cluster inside "F. nucleatum" is `divergent` and stays under
               nucleatum (Ammal, 26 Oct 2026): Bracken can only call those reads F. nucleatum.

groups(...)    Labels whose anchors are joined at >= group_ani (single linkage over anchor genome
               pairs) form one ANI group, named "<lowest-taxid member>_group". 92% (Ammal, 26 Oct
               2026) puts F. nucleatum, animalis, polymorphum, vincentii, hwasookii and canifelinum
               in one group; 93% would leave polymorphum out. The scorer uses the group only as a
               fallback level for cells that species-level reads cannot resolve.
"""
import statistics
from collections import defaultdict

UNNAMED = "_sp._"


def genus_of(label):
    return label.split("_", 1)[0]


def _clusters(members, ani, threshold):
    parent = {m: m for m in members}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, a in enumerate(members):
        for b in members[i + 1:]:
            if ani.get((a, b), 0.0) >= threshold:
                parent[find(a)] = find(b)
    comps = defaultdict(list)
    for m in members:
        comps[find(m)].append(m)
    return list(comps.values())


def anchors(genomes, ani, anchor_ani=95.0):
    """genomes: ordered [(accession, label)], first genome of a label = its seed. -> {label: set(acc)}"""
    members = defaultdict(list)
    for acc, label in genomes:
        members[label].append(acc)
    out = {}
    for label, accs in members.items():
        seed = accs[0]
        out[label] = set(next(c for c in _clusters(accs, ani, anchor_ani) if seed in c))
    return out


def relabel(genomes, ani, anchor_ani=95.0, margin=1.0):
    """-> list of dicts (accession, filed_label, label, action, score_own, best_label, score_best);
    label is '' for a dropped genome."""
    anc = anchors(genomes, ani, anchor_ani)
    members = defaultdict(list)
    for acc, label in genomes:
        members[label].append(acc)
    out = []
    for acc, L in genomes:
        scores = {}
        for M, A in anc.items():
            if M != L and UNNAMED in M:
                continue
            v = [ani.get((acc, h), 0.0) for h in A if h != acc]
            if v:
                scores[M] = statistics.median(v)
        own = scores.get(L, 0.0)
        best = max(sorted(scores), key=lambda M: scores[M]) if scores else L
        if acc == members[L][0] and UNNAMED not in L:
            action, new = "keep_seed", L
        elif scores.get(best, 0.0) < anchor_ani:
            action, new = ("keep_singleton" if len(members[L]) == 1 else "divergent"), L
        elif best == L:
            action, new = "keep", L
        elif scores[best] - own >= margin:
            action, new = "relabel", best
        else:
            action, new = "drop_ambig", ""
        out.append(dict(accession=acc, filed_label=L, label=new, action=action,
                        score_own=own, best_label=best, score_best=scores.get(best, 0.0)))
    return out


def groups(genomes, ani, taxid_of, anchor_ani=95.0, group_ani=92.0):
    """genomes AFTER relabelling, ordered. -> {label: group_name} for labels in a group of >= 2."""
    anc = anchors(genomes, ani, anchor_ani)
    labels = sorted(anc)
    parent = {L: L for L in labels}

    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x

    for i, L in enumerate(labels):
        for M in labels[i + 1:]:
            if any(ani.get((a, b), 0.0) >= group_ani for a in anc[L] for b in anc[M]):
                parent[find(L)] = find(M)
    comps = defaultdict(list)
    for L in labels:
        comps[find(L)].append(L)
    out = {}
    for comp in comps.values():
        if len(comp) < 2:
            continue
        def key(L):
            t = taxid_of.get(L, "")
            return (0, int(t), L) if t.isdigit() else (1, 0, L)
        name = min(comp, key=key) + "_group"
        for L in comp:
            out[L] = name
    return out
