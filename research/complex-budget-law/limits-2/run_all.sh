#!/bin/bash
# Replay every check in limits-2 from its own directory and report the PASS_ marker of each.
# Usage: bash run_all.sh [--quick]     (--quick skips the runs longer than ~20 s)
# Python: system python3 for stdlib scripts; uv (Python 3.12) for scripts that need numpy/scipy/networkx/ortools.
# One worker: BLAS/OpenMP threads are pinned to 1 and CP-SAT runs with num_workers = 1.
# Some inputs live elsewhere in this branch: theta-partial/xr imports ../../limits-1/complex-dual/{xk,xg} (expXK, expXG and
# their graph files); contract-slack/xv and complex-search/xp read ../../limits-1/complex-dual/xg/graphs_small.json.
# pending/ holds a placeholder only and has no checks.
set -u
cd "$(dirname "$0")"
QUICK=${1:-}
UV="uv run -q --no-project --python 3.12"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
fail=0; n=0
run() {   # run <dir> <slow:0|1> <cmd...>
    local dir=$1 slow=$2; shift 2
    if [ "$QUICK" = "--quick" ] && [ "$slow" = 1 ]; then echo "SKIP  $dir: $*"; return; fi
    local t0=$SECONDS out
    out=$(cd "$dir" && "$@" 2>&1 | tail -1); n=$((n + 1))
    if [[ "$out" == PASS* ]]; then echo "ok    $dir: $* -> $out ($((SECONDS - t0)) s)"
    else echo "FAIL  $dir: $* -> $out"; fail=1; fi
}
P="python3 -I -B"
NP="$UV --with numpy python -I -B"
G=../../../limits-1/complex-dual/xg/graphs_small.json
# theta-partial (expXR)
run theta-partial/xr 1 $UV --with numpy --with scipy python -I -B check_xr.py
run theta-partial/xr 1 $UV --with numpy --with scipy python -I -B check_xr2.py
# exact-instances (expXQ, expXU)
run exact-instances/xq 0 $P check_json.py schedule9.json
run exact-instances/xq 0 $NP verify.py 3 cyc9.txt dense
run exact-instances/xq 0 env PERTURB=one $NP verify.py 3 cyc9.txt dense
run exact-instances/xq 0 bash -c "$NP mkgraph.py 3 && $NP cyc.py 9"
run exact-instances/xq 0 bash -c "$NP mkgraph.py 5 && $NP cyc.py 25"
run exact-instances/xq 0 $NP verify.py 5 cyc25.txt
run exact-instances/xq 1 bash -c "$NP mkgraph.py 7 && $NP cyc.py 49; rc=\$?; rm -f graph49.txt; exit \$rc"
run exact-instances/xu 0 $P summarize.py
# contract-slack (expL8, expXW, expXV)
run contract-slack/l8 0 $NP check_l8.py
run contract-slack/xw 0 $P xw.py 60
run contract-slack/xv 0 $UV --with numpy --with scipy --with networkx python -I -B flow_check.py
# complex-search (expXP)
run complex-search/xp 0 $NP spectral_boot.py
run complex-search/xp 0 $P ladder_summary.py
run complex-search/xp 0 env PASS_F=1e-10 $NP cmds.py $G PG3 8 0
run complex-search/xp 1 env PASS_F=1e-6 $NP cmds.py $G PG5 16 4 4000
# full-pg16 (expXX)
run full-pg16/xx 0 $UV --with numpy --with scipy python -I -B markov17.py
run full-pg16/xx 0 bash -c "$UV --with numpy --with ortools python -I -B sim16.py seed_11_91.txt sched909.json; rc=\$?; rm -f code16.pkl; exit \$rc"
echo "checks run: $n, exit $fail"
exit $fail
