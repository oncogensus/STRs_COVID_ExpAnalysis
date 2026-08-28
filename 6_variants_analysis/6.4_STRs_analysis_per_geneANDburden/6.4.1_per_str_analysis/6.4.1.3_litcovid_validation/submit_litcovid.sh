#!/usr/bin/env bash
# submit_litcovid.sh
cd "$(dirname "${BASH_SOURCE[0]}")"
qsub litcovid.pbs
