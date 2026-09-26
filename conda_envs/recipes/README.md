# Local conda recipes for the batch-correction environment

ConQuR (the ivartb/ConQuR_par fork) and cqrReg 1.2.1 are on no conda channel the pipeline uses
(checked 2026-09-26: absent from conda-forge and bioconda; cqrReg exists only on Anaconda's `r`
channel, built against Anaconda's R). cqrReg 1.2.1 is the current CRAN release. These recipes build them as conda packages from the same
checksum-verified sources that scripts/install_conqur.R used, against the exact r-base, compiler
and dependency builds of conda_envs/lock/batch_correction_env.linux-64.txt, so the compiled code
matches the environment the pipeline was validated with.

The built packages are committed under conda_envs/local_channel/ and listed in the batch lockfile;
Nextflow then builds the batch environment itself, like every other environment.
Rebuild with scripts/build_local_conda_packages.sh (Linux x86_64, network for conda-forge).
