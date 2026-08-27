#!/usr/bin/env bash
# submit_dbscan_subset.sh
cd "$(dirname "${BASH_SOURCE[0]}")"
qsub dbscan_subset.pbs
