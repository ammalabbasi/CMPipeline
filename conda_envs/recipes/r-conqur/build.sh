#!/bin/bash
# The same Remote* fields scripts/install_conqur.R appended for the offline tarball install.
set -euo pipefail
# Keep equal to params.conqur_sha (main.nf) and the sha in meta.yaml (tests/test_local_conda_channel.py).
CONQUR_SHA=ff233085cedc36a24382038067a4cece3e94e0a0
cat >> DESCRIPTION <<DESC
RemoteType: github
RemoteHost: api.github.com
RemoteUsername: ivartb
RemoteRepo: ConQuR_par
RemoteRef: ${CONQUR_SHA}
RemoteSha: ${CONQUR_SHA}
GithubRepo: ConQuR_par
GithubUsername: ivartb
GithubSHA1: ${CONQUR_SHA}
DESC
${R} CMD INSTALL --build .
