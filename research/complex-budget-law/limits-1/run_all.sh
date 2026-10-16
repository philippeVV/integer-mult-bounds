#!/bin/bash
# Replay every check in limits-1 from its own directory and report the PASS_ marker of each.
# Usage: bash run_all.sh [--quick]     (--quick skips the runs longer than ~20 s)
# Python: system python3 for stdlib scripts; uv (Python 3.12) for scripts that need numpy/scipy/numba/ortools/networkx/
# cvxpy/clarabel. One worker: BLAS/OpenMP/numba threads are pinned to 1 and CP-SAT runs with num_workers = 1.
# Some inputs live elsewhere in this branch: lattices/xa/price.py imports ../../../framework-ceiling and closures-3 copies.
set -u
cd "$(dirname "$0")"
QUICK=${1:-}
UV="uv run -q --no-project --python 3.12"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
NUMBA_CACHE_DIR=$(mktemp -d); export NUMBA_CACHE_DIR; trap 'rm -r "$NUMBA_CACHE_DIR"' EXIT
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
# single-letter (expXB)
run single-letter 0 $P xb.py
# complex-dual (expXG, expXK)
run complex-dual/xg 0 $UV --with numpy --with scipy python -I -B check_xg.py
run complex-dual/xk 0 $UV --with numpy --with scipy python -I -B check_xk.py
run complex-dual/xk 0 $UV --with numpy --with scipy python -I -B sic_schur.py
# bootstrap (expXE; expXL pending)
run bootstrap/xe 0 $P check_boot.py
run bootstrap/xe 0 env GF=graphs_small.json $UV --with networkx python -I -B boot3.py PG3 100 8
for d in 4 12; do run bootstrap/xe 0 $P bootrand.py $d 1000 10; done
# suppliers (expXF, expXJ; expXN pending)
run suppliers/xf 0 $P xf.py
run suppliers/xj 0 $P psi_check.py
run suppliers/xj 0 $P control.py
run suppliers/xj 0 $UV --with numpy python -I -B odd_singer.py 3 dense
run suppliers/xj 0 $UV --with numpy python -I -B odd_singer.py 5
run suppliers/xj 1 $UV --with numpy python -I -B odd_singer.py 7
# lattices (expXA, expXH)
run lattices/xa 0 $UV --with numpy python -I -B lor.py
run lattices/xa 0 $UV --with numpy python -I -B scale.py
run lattices/xa 0 $UV --with numpy python -I -B rankW.py
run lattices/xa 0 $UV --with numpy python -I -B bw16.py
run lattices/xa 0 $UV --with numpy python -I -B price.py
run lattices/xa 0 $UV --with numpy --with cvxpy --with clarabel python -I -B theta_lzd.py 9 11
run lattices/xh 0 $UV --with numpy --with ortools python -I -B part1.py 12 3 e8pad 60
run lattices/xh 0 $UV --with numpy --with ortools python -I -B part1.py 12 3 rand 10 5
run lattices/xh 0 $UV --with numpy --with ortools python -I -B part2.py
for k in 3 4 5 6 7; do run lattices/xh 0 $UV --with numpy --with ortools python -I -B part2b.py gauss $k 4 k; done
for k in 5 7; do run lattices/xh 0 $UV --with numpy --with ortools python -I -B part2b.py eis $k 4 k; done
run lattices/xh 0 $UV --with numpy python -I -B price2.py
# seeds-negative (expXD, expXM, expRJ; expRB/expRF are in ../closures-3/complex-coding-seeds)
run seeds-negative/xd 0 $UV --with numpy --with numba python -I -B check_xd.py
run seeds-negative/xd 1 $UV --with numpy --with numba python -I -B h5.py 5
run seeds-negative/xm 0 $UV --with numpy python -I -B fam.py 0
run seeds-negative/xm 0 $UV --with numpy python -I -B fam.py 1
run seeds-negative/xm 0 $P agree.py
run seeds-negative/xm 0 $UV --with numpy --with scipy python -I -B lm.py m1:120 3
run seeds-negative/rj 0 $P summarize.py
echo "checks run: $n, exit $fail"
exit $fail
