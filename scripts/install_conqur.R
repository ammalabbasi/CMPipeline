#!/usr/bin/env Rscript
# Installs the two batch-correction packages conda cannot provide, at EXACT versions, into the
# active environment (built from conda_envs/batch_correction_env.yml; see
# conda_envs/build_envs.sbatch, which runs this for you):
#
#   * cqrReg 1.2.1        -- a ConQuR import that is not on conda-forge (from the CRAN archive)
#   * ConQuR              -- the ivartb/ConQuR_par fork, which is what
#                            scripts/2500703_batch_correction_normalization.r was written against
#                            (it fixes a "batchid not found" error in ConQuR's foreach workers)
#
# Everything else ConQuR imports comes from conda, so both installs use dependencies = FALSE:
# nothing is fetched from CRAN at whatever version it serves today (audit D03, 2026-09-23).
# The earlier default, wdl2459/ConQuR@master, named a branch that does not exist upstream.

repo  <- Sys.getenv("CONQUR_REPO", unset = "ivartb/ConQuR_par")
ref   <- Sys.getenv("CONQUR_REF",  unset = "ff233085cedc36a24382038067a4cece3e94e0a0")
force <- nzchar(Sys.getenv("CONQUR_FORCE"))
cqr_version <- "1.2.1"
# Offline source tarballs (audit D03): used first when present, checksum-verified against SHA256SUMS,
# so the env can be rebuilt without CRAN or GitHub (cqrReg 1.2.1 leaves CRAN's main index as soon as
# a newer version appears). Falls back to the network when the directory is absent.
offline <- Sys.getenv("CONQUR_OFFLINE_DIR",
                      unset = "/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/packages/conqur_offline")
offline_file <- function(name) {
    path <- file.path(offline, name)
    sums <- file.path(offline, "SHA256SUMS")
    if (!file.exists(path) || !file.exists(sums)) return(NULL)
    want <- sub(" .*", "", grep(paste0(" ", name, "$"), readLines(sums), value = TRUE))
    got <- sub(" .*", "", system2("sha256sum", path, stdout = TRUE))
    if (length(want) != 1 || !identical(want, got)) stop("checksum mismatch for ", path)
    path
}

if (!requireNamespace("remotes", quietly = TRUE))
    stop("r-remotes is missing: build the env from conda_envs/batch_correction_env.yml first")

if (!requireNamespace("cqrReg", quietly = TRUE) ||
        as.character(utils::packageVersion("cqrReg")) != cqr_version) {
    tb <- offline_file(paste0("cqrReg_", cqr_version, ".tar.gz"))
    if (!is.null(tb)) {
        message("installing cqrReg from ", tb)
        install.packages(tb, repos = NULL, type = "source")
    } else {
        remotes::install_version("cqrReg", version = cqr_version, dependencies = FALSE,
                                 upgrade = "never", repos = "https://cloud.r-project.org")
    }
}

have_pin <- requireNamespace("ConQuR", quietly = TRUE) &&
    identical(utils::packageDescription("ConQuR")$RemoteSha, ref)
if (have_pin && !force) {
    message("ConQuR already installed at ", repo, "@", ref)
} else {
    tb <- offline_file(paste0(basename(repo), "-", ref, ".tar.gz"))
    if (!is.null(tb)) {
        # GitHub's codeload tarball of exactly this commit (checksum verified). Record the same
        # Remote* fields install_github writes, so the module's RemoteSha check sees the pin.
        message("installing ConQuR from ", tb)
        src <- tempfile(); dir.create(src); utils::untar(tb, exdir = src)
        pkg <- list.dirs(src, recursive = FALSE)[1]
        write(c("RemoteType: github", "RemoteHost: api.github.com",
                paste0("RemoteUsername: ", dirname(repo)), paste0("RemoteRepo: ", basename(repo)),
                paste0("RemoteRef: ", ref), paste0("RemoteSha: ", ref),
                paste0("GithubRepo: ", basename(repo)), paste0("GithubUsername: ", dirname(repo)),
                paste0("GithubSHA1: ", ref)),
              file = file.path(pkg, "DESCRIPTION"), append = TRUE)
        install.packages(pkg, repos = NULL, type = "source")
    } else {
        remotes::install_github(paste0(repo, "@", ref), dependencies = FALSE, upgrade = "never", force = TRUE)
    }
}

# Verify: the pinned commit, and every ConQuR import loadable.
stopifnot(identical(utils::packageDescription("ConQuR")$RemoteSha, ref))
imports <- c("quantreg", "cqrReg", "glmnet", "dplyr", "doParallel", "gplots", "vegan", "ade4",
             "compositions", "randomForest", "ROCR", "ape", "GUniFrac", "fastDummies")
missing <- imports[!vapply(imports, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("ConQuR imports missing from the env: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(ConQuR))
message("ConQuR ", as.character(utils::packageVersion("ConQuR")), " ready (", repo, "@", ref,
        "); cqrReg ", as.character(utils::packageVersion("cqrReg")))
