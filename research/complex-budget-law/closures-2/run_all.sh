#!/bin/bash
# Replay every check in closures-2 from its own lane directory and report the PASS_ marker of each.
# Usage: bash run_all.sh [--quick]     (--quick skips the runs longer than ~1 min)
# Python: system python3 for stdlib scripts; uv (Python 3.12) for scripts that need numpy/ortools/python-sat/python-flint.
set -u
cd "$(dirname "$0")"
QUICK=${1:-}
UV="uv run -q --no-project --python 3.12"
fail=0
run() {   # run <lane> <slow:0|1> <cmd...>
    local lane=$1 slow=$2; shift 2
    if [ "$QUICK" = "--quick" ] && [ "$slow" = 1 ]; then echo "SKIP  $lane: $*"; return; fi
    local t0=$SECONDS out
    out=$(cd "$lane" && "$@" 2>&1 | tail -1)
    if [[ "$out" == PASS* ]]; then echo "ok    $lane: $* -> $out ($((SECONDS - t0)) s)"
    else echo "FAIL  $lane: $* -> $out"; fail=1; fi
}
P="python3 -I -B"
run flat-tower 0 $P chain.py
run flat-tower 0 $P levels.py
run flat-tower 0 $P stage.py
run flat-tower 1 $P price_flat.py
run flat-tower 0 $P price_n.py
run zll-recursion 0 $P check_qg.py
run zll-recursion 0 $P zll_scale.py
run zll-recursion 0 $P check_qm.py
run zll-recursion 1 $P price_qm.py
run sibling-helpers 0 $P census.py
run sibling-helpers 1 $P readercut.py
run sibling-helpers 1 $P price.py
run p120-rank8-rigidity 0 $P rupcheck.py
run p120-rank8-rigidity 0 $UV --with numpy --with python-sat python -I -B nondeg_cnf.py
run p120-rank8-rigidity 0 $UV --with numpy --with python-flint python -I -B tangent.py
run p120-rank8-rigidity 0 $UV --with numpy --with ortools python -I -B alpha.py
run p120-rank8-rigidity 0 $UV --with numpy --with ortools python -I -B unions.py 0 3
run p120-rank8-rigidity 0 $P floor.py
run p120-rank8-rigidity 0 $P halfint.py
run p120-rank8-rigidity 0 $P score.py
run p120-rank8-rigidity 0 $P price_rank.py
run tm-pass-programs 0 $P lagbound.py
run tm-pass-programs 1 $UV --with python-sat python -I -B replay_small.py
for fam in "2 3 dot" "2 4 herm" "2 4 dot" "2 4 dot plus1" "3 3 dot plus1"; do
    run field-line-large-ports 0 $UV --with numpy --with ortools python -I -B exact.py $fam
done
for fam in "4 3 herm" "4 3 dot" "3 4 dot"; do
    run field-line-large-ports 1 $UV --with numpy --with ortools python -I -B exact.py $fam
done
run field-line-large-ports 0 $P kappa.py
run omega-h-floor-false 0 $P verify.py
run small-units-h5 0 $P score.py
run small-units-h5 0 $P verify_hit.py
exit $fail
