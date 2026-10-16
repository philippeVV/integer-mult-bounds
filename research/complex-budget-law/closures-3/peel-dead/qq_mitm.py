#!/usr/bin/env python3
"""expQQ: meet-in-the-middle enumeration of ALL all-role permutations rho realisable on a fixed two-port graph
(upstream stronger-rank-screens.md section 4 method: prefix A, inverse suffix B^-1, BA = P iff the column multisets
agree), then for every NON-routable rho (the only ones not already excluded by char-0 routing):
  CP-SAT  min #changed edges, scalar frames, any climbs      (Theorem A at m = 1 says = W)
  CP-SAT  same, with at most one climb per role route         (single climb, scalar)
usage: qq_mitm.py fano10 | rand W n seed [count]
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys
import time
import random
from itertools import product
from qq_core import OPS, routable, routes

FANO10 = ((0, 3), (1, 4), (2, 5), (0, 1), (4, 2), (0, 6), (4, 7), (0, 4), (5, 6), (0, 8), (5, 9), (0, 5), (7, 9), (3, 8))


def apply_cols(cols, a, b, f):
    out = []
    for c in cols:
        xa, xb = c >> a & 1, c >> b & 1
        ya, yb = f(xa, xb)
        c = c & ~((1 << a) | (1 << b)) | (ya << a) | (yb << b)
        out.append(c)
    return tuple(out)


def half(W, gates):
    states = {tuple(1 << j for j in range(W))}
    for a, b in gates:
        states = {apply_cols(s, a, b, f) for s in states for f in OPS}
    return states


def all_perms(W, pairs):
    k = len(pairs) // 2
    A = half(W, pairs[:k])
    Binv = half(W, tuple(reversed(pairs[k:])))   # B^-1 = G_k^-1 ... G_n^-1 (ops closed under inverse)
    byA = {}
    for s in A:
        key = tuple(sorted(s))
        byA.setdefault(key, set()).add(tuple(key.index(c) for c in s))   # sigma_A(j)
    rhos = set()
    for s in Binv:
        key = tuple(sorted(s))
        if key not in byA:
            continue
        pos = {c: i for i, c in enumerate(s)}                            # column value -> index in B^-1
        for sA in byA[key]:
            rhos.add(tuple(pos[key[sA[j]]] for j in range(W)))
    return rhos, len(A), len(Binv)


def check_rho_scalar(W, pairs):
    from qq_core import scalar_min_changes
    return scalar_min_changes


def cpsat_min(W, pairs, rho, single_climb):
    from ortools.sat.python import cp_model
    n = len(pairs)
    R = routes(W, pairs)
    N = 2 * W + n
    K = N + 2
    md = cp_model.CpModel()
    p = [md.NewIntVar(-K, K, f'p{v}') for v in range(N)]
    bs = []
    for r in R:
        rb = []
        for i in range(len(r) - 1):
            b = md.NewBoolVar('')
            md.Add(p[r[i]] == p[r[i + 1]]).OnlyEnforceIf(b.Not())
            rb.append(b)
        if single_climb:
            md.Add(sum(rb) <= 1)
        bs += rb
    for w in range(W):
        md.Add(p[W + n + rho[w]] == p[w] + 1)
    md.Minimize(sum(bs))
    sv = cp_model.CpSolver()
    sv.parameters.num_workers = 2
    sv.parameters.max_time_in_seconds = 60
    st = sv.Solve(md)
    if st == cp_model.INFEASIBLE:
        return 'infeasible'
    assert st == cp_model.OPTIMAL, sv.StatusName(st)
    return int(round(sv.ObjectiveValue()))


def study(W, pairs, label, maxcheck=200):
    t0 = time.time()
    rhos, nA, nB = all_perms(W, pairs)
    nonr = [r for r in sorted(rhos) if not routable(W, pairs, r)]
    print(f'{label}: W={W} n={len(pairs)} halves {nA}/{nB}; realisable rho {len(rhos)}; non-routable {len(nonr)}'
          f' [{time.time()-t0:.1f}s]', flush=True)
    worst_any = worst_one = None
    for r in nonr[:maxcheck]:
        a = cpsat_min(W, pairs, r, False)
        o = cpsat_min(W, pairs, r, True)
        if a < W or (o != 'infeasible' and o < W):
            print('   VIOLATION', r, a, o, flush=True)
        worst_any = a if worst_any is None else min(worst_any, a)
        if o != 'infeasible':
            worst_one = o if worst_one is None else min(worst_one, o)
    if nonr:
        print(f'   checked {min(len(nonr), maxcheck)} non-routable rho: min changed edges (any climbs) = {worst_any},'
              f' (single climb) = {worst_one}  (W = {W})  [{time.time()-t0:.1f}s]', flush=True)
    return rhos, nonr


def main():
    if sys.argv[1] == 'fano10':
        rhos, nonr = study(10, FANO10, 'upstream sec.4 balanced Fano graph')
        # upstream terminal obligations: sources 0,1,2 -> 7,0,3
        hit = sum(1 for r in rhos if r[:3] == (7, 0, 3))
        print('   rho with rho(0),rho(1),rho(2) = 7,0,3:', hit)
        print('PASS_EXPQQ_MITM' if (len(rhos), len(nonr), hit) == (14064, 0, 0) else 'FAIL_EXPQQ_MITM')
    else:
        W, n, seed = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
        count = int(sys.argv[5]) if len(sys.argv) > 5 else 1
        rng = random.Random(seed)
        tot_nonr = 0
        for t in range(count):
            pairs = tuple(tuple(rng.sample(range(W), 2)) for _ in range(n))
            _, nonr = study(W, pairs, f'rand seed={seed} #{t}', maxcheck=40)
            tot_nonr += len(nonr)
        print('total non-routable rho:', tot_nonr)


if __name__ == '__main__':
    main()
