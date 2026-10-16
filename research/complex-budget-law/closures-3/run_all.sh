#!/bin/bash
# Replay every check in closures-3 from its own directory and report the PASS_ marker of each.
# Usage: bash run_all.sh [--quick]     (--quick skips the runs longer than ~1 min)
# Python: system python3 for stdlib scripts; uv (Python 3.12) for scripts that need numpy/scipy/ortools/networkx/mpmath/
# cvxpy/python-flint. C programs are compiled by the lane's own shell driver (gcc -O2) into a temporary directory.
set -u
cd "$(dirname "$0")"
QUICK=${1:-}
UV="uv run -q --no-project --python 3.12"
fail=0
run() {   # run <dir> <slow:0|1> <cmd...>
    local dir=$1 slow=$2; shift 2
    if [ "$QUICK" = "--quick" ] && [ "$slow" = 1 ]; then echo "SKIP  $dir: $*"; return; fi
    local t0=$SECONDS out
    out=$(cd "$dir" && "$@" 2>&1 | tail -1)
    if [[ "$out" == PASS* ]]; then echo "ok    $dir: $* -> $out ($((SECONDS - t0)) s)"
    else echo "FAIL  $dir: $* -> $out"; fail=1; fi
}
P="python3 -I -B"
# peel-dead
run peel-dead 0 $P hyp_drop.py
run peel-dead 0 $UV --with mpmath python -I -B lemma_p.py
run peel-dead 0 $UV --with numpy --with scipy python -I -B qq_enum.py 3 5
run peel-dead 1 $P qq_mitm.py fano10
run peel-dead 0 $UV --with ortools python -I -B rect_lemma.py
run peel-dead 0 $UV --with ortools python -I -B rect_max.py
run peel-dead 0 $P conic_search.py 2 2 500 1 2 rank1 12 4
run peel-dead 0 $P qy_check.py all
for L in 2,2,2 1,2,2,1 2,3,2; do run peel-dead 0 $UV --with ortools python -I -B hyper_cg.py $L 1000; done
run peel-dead 0 $UV --with ortools python -I -B ctl_realise.py 2,2,2 '2-3;3-2' 192
run peel-dead 0 $UV --with ortools python -I -B route_lek.py 4 2,2,2
run peel-dead 0 $UV --with ortools python -I -B classgraph_sat.py --claim far
run peel-dead 0 $UV --with ortools python -I -B classgraph_sat.py --claim far --cuts none --expect FEASIBLE
run peel-dead 0 $UV --with numpy python -I -B rh_brute_check.py 3 2 4 2 2 3
# large-ports-trace-bound
run large-ports-trace-bound 0 $P trace_bound.py
run large-ports-trace-bound 0 $UV --with numpy --with cvxpy python -I -B theta_local.py
run large-ports-trace-bound 0 $UV --with numpy python -I -B certify.py
run large-ports-trace-bound 0 $UV --with numpy --with networkx python -I -B obstruct.py
run large-ports-trace-bound 0 $P price.py
run large-ports-trace-bound/qa 0 $UV --with numpy python -I -B verify.py
# laminar-decoders
run laminar-decoders 0 $P dfs.py
run laminar-decoders 0 $P minbow.py
run laminar-decoders 0 $UV --with ortools python -I -B lam.py
run laminar-decoders 0 $P maxlam.py
run laminar-decoders 0 $UV --with numpy --with scipy python -I -B span_milp.py
run laminar-decoders 0 $P price.py
# endpoint-labels-both-banks
run endpoint-labels-both-banks 0 $UV --with numpy python -I -B eqv_motif_y.py 12
run endpoint-labels-both-banks 0 $UV --with numpy python -I -B eqv_motif_y.py 40 neg
run endpoint-labels-both-banks 0 $UV --with networkx python -I -B match_p112.py
run endpoint-labels-both-banks 0 $P price.py
run endpoint-labels-both-banks 0 $UV --with networkx python -I -B flat.py
run endpoint-labels-both-banks 0 $UV --with networkx python -I -B replay_y.py
run endpoint-labels-both-banks 0 bash bfs_small.sh
# char2-tanner-family (family record)
run char2-tanner-family/qd 0 $P verify_cut.py
run char2-tanner-family/qe 0 $UV --with numpy python -I -B hilb.py
run char2-tanner-family/qe 0 $UV --with numpy python -I -B verify_sched.py
run char2-tanner-family/qe 0 $UV --with numpy python -I -B colrep_bound.py
run char2-tanner-family/qu 0 $UV --with numpy --with ortools python -I -B divcheck.py
run char2-tanner-family/qk 0 $UV --with numpy python -I -B qk_singer.py 2 1 0
run char2-tanner-family/qk 0 $UV --with numpy python -I -B qk_singer.py 4 2 0
run char2-tanner-family/qk 0 $UV --with numpy python -I -B qk_singer.py 6 3 0
run char2-tanner-family/qi 0 $UV --with numpy --with scipy python -I -B flow.py 3 1 0 8 0
run char2-tanner-family/qv 0 $P table.py
run char2-tanner-family/ny 1 $UV --with numpy python -I -B check_telescope.py
run char2-tanner-family/qr 0 $UV --with numpy python -I -B check_lemma.py
run char2-tanner-family/qr 0 $UV --with numpy --with python-flint python -I -B families.py 0 same 2
run char2-tanner-family/qr 0 $UV --with numpy --with python-flint python -I -B ray.py same 2 F1_Droot_r1
run char2-tanner-family/qr 0 $UV --with numpy --with python-flint python -I -B cutrank.py
run char2-tanner-family/ra 0 $P check_flow.py q16.json cert_flow_q16.json.gz
run char2-tanner-family/ra 0 $P check_flow.py s784.json cert_flow_s784.json.gz
run char2-tanner-family/l6 0 $P estimates.py
# complex-coding-seeds
run complex-coding-seeds/rb 0 $UV --with numpy --with galois python -I -B c1.py 3
run complex-coding-seeds/rb 0 $P table.py
run complex-coding-seeds/rf 0 $UV --with numpy --with scipy python -I -B markov.py
run complex-coding-seeds/rf 0 $UV --with numpy --with galois python -I -B equiv.py 8 3000 1
run complex-coding-seeds/rf 0 $P table.py
# small-units-lp
run small-units-lp 0 $UV --with python-flint python -I -B spread.py
run small-units-lp 1 $UV --with numpy --with python-flint python -I -B idem77.py
run small-units-lp 0 $UV --with numpy --with python-flint python -I -B rpub.py 5 7 4
# decoder-law-subsets
run decoder-law-subsets/qj 0 bash ladder.sh
run decoder-law-subsets/qj 0 $UV --with numpy python -I -B check.py
run decoder-law-subsets/r3 0 $UV --with numpy python -I -B csys_check.py 1
run decoder-law-subsets/r3 1 $UV --with numpy python -I -B members.py 1 9
run decoder-law-subsets/nz 0 $UV --with ortools python -I -B cluster.py 900
run decoder-law-subsets/nz 1 $P certhc.py 1
exit $fail
