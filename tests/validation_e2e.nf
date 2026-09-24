// End-to-end test of the alignment-validation processes on SIMULATED reads from public genomes.
// Not part of the pipeline. Runs the real processes (species list -> reference build with
// storeDir cache -> batched streamed alignment -> scoring -> mask -> level comparison) and
// writes their outputs to --outdir for tests/check_validation_e2e.py to check against the truth.
//
//   nextflow run tests/validation_e2e.nf -c tests/validation_e2e.config \
//       --e2e_dir <dir with samples.csv, species_table.tsv> --outdir <out>
nextflow.enable.dsl=2

params.e2e_dir = null
params.outdir = null
params.scripts = "${projectDir}/../scripts"
params.validation_env = "${projectDir}/../conda_envs/validation_env.yml"
params.validation_dir = "${params.outdir}/06_VALIDATION"
params.validation_cache_dir = "${params.outdir}/validation_cache"
params.validation_groups = "${projectDir}/../ref/validation_groups.tsv"
params.validation_taxdb = null
params.validation_min_support = 25
params.validation_min_bins_ratio = 0.5
params.validation_min_identity = 0.98
params.validation_min_unique_frac = 0.02
params.validation_max_genomes = 5
params.validation_dedup_ani = 99.5
params.validation_max_missing_frac = 0.5
params.validation_batch_size = 2
params.validation_export_reads = true

// Same dependency digests main.nf passes (audit C03), so this harness also tests -resume behaviour.
def dig(List names) { names.collect { file("${params.scripts}/${it}").text.md5() }.join(";") }

include { validationSpeciesList; buildValidationRef; alignForValidation; applyValidationMask; compareDecontamLevels } from '../Modules/alignment_validation.nf'

workflow {
    def table = file("${params.e2e_dir}/species_table.tsv", checkIfExists: true)
    def groups = file(params.validation_groups, checkIfExists: true)
    Channel.fromPath("${params.e2e_dir}/samples.csv").splitCsv(header: true)
        .map { r -> tuple(r.sample, file(r.reads, checkIfExists: true)) }
        .collate(params.validation_batch_size)
        .map { batch -> tuple(batch.collect { it[0] }, batch.collect { it[1] }) }
        .set { BATCHES }
    validationSpeciesList(table, file(params.validation_taxdb, checkIfExists: true), dig(["validation_species_list.py", "validation_common.py", "build_validation_ref.py"]))
    buildValidationRef(validationSpeciesList.out.species)
    alignForValidation(BATCHES, buildValidationRef.out.ref.first(), table, groups, dig(["validation_score.py", "validation_common.py"]))
    applyValidationMask(table, alignForValidation.out.stats.flatten().collect(), groups, dig(["apply_validation_mask.py", "validation_common.py"]))
    compareDecontamLevels(file("${params.e2e_dir}/species.contaminants_final.tsv"),
                          file("${params.e2e_dir}/genus.contaminants_final.tsv"), dig(["compare_decontam_levels.py", "validation_common.py"]))
}
