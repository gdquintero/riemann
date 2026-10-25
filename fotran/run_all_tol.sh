#!/usr/bin/env bash
# Equivalent of ../run_all_tol.py: the three tolerance experiments, one after the other
# (the computation is memory bound, so running them at the same time does not help).
set -euo pipefail
cd "$(dirname "$0")"
make -s
for tol in 0.01 0.001 0.0001; do
    ./riemann_tol "riemann_${tol}.nml"
done
