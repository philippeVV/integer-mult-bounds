#!/usr/bin/env python3
"""expQQ exhaustive small-gadget search (Q-peel stop test, small sizes).

For every two-port gate sequence on W roles with n gates (up to role relabelling) and every all-role permutation
rho realisable by GL(2,F2) gates on it:
  T  scalar frames, any climbs: min #changed physical edges (must be >= W by Theorem A at m = 1);
  M  two frame classes per role route, MONOTONE (single climb, A..A B..B): exact gluing by union-find,
     lower bound s/m >= hard + LP (fractional packing of unbalanced cycles), upper bound hard + scalar;
  X  two frame classes per role route, ANY pattern (A B A ... allowed, cost = #changes x rank): same bounds.
Prints the minimum of s/(W m) found (LP lower bound and scalar upper bound) per (W, n).
usage: qq_enum.py WMAX NMAX [XMAX_PATTERNS]
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys
import time
from itertools import permutations, product, combinations
from qq_core import (achievable_perms, routes, routable, scalar_min_changes, route_patterns, solve_pattern,
                     soft_scalar_min, soft_lp)


def canon_seqs(W, n):
    pairs = list(combinations(range(W), 2))
    seen = set()
    perms = list(permutations(range(W)))
    for seq in product(pairs, repeat=n):
        best = min(tuple(tuple(sorted((p[a], p[b]))) for a, b in seq) for p in perms)
        if best in seen:
            continue
        seen.add(best)
        yield best


def pattern_bounds(W, n, R, rho, monotone, cap):
    pats = [list(route_patterns(len(R[w]) - 2, monotone)) for w in range(W)]
    total = 1
    for p in pats:
        total *= len(p)
    if total > cap:
        return None
    lo = hi = float('inf')
    feas = 0
    for pat in product(*pats):
        r = solve_pattern(W, n, R, rho, pat)
        if r is None:
            continue
        feas += 1
        hard, soft = r
        if hard >= W:
            lo = min(lo, hard); hi = min(hi, hard + sum(s[3] for s in soft)); continue
        if not soft:
            lo = min(lo, hard); hi = min(hi, hard); continue
        sc = hard + soft_scalar_min(soft)
        hi = min(hi, sc)
        if sc < W:
            lo = min(lo, sc)  # scalar solution is a real frame solution (times I_m)
            print('   HIT scalar', W, n, rho, pat, hard, soft, flush=True)
            continue
        lp = hard + soft_lp(soft)
        lo = min(lo, lp)
        if lp < W - 1e-9:
            print('   FLAG lp<W', W, n, rho, pat, hard, soft, round(lp, 4), flush=True)
    return lo, hi, feas, total


def main():
    WMAX, NMAX = int(sys.argv[1]), int(sys.argv[2])
    XCAP = int(sys.argv[3]) if len(sys.argv) > 3 else 1 << 13
    ok = True
    for W in range(2, WMAX + 1):
        for n in range(1, NMAX + 1):
            t0 = time.time()
            nseq = ninst = nnr = 0
            tmin = W
            Mlo = Mhi = Xlo = Xhi = float('inf'); Xskip = 0
            for seq in canon_seqs(W, n):
                nseq += 1
                R = routes(W, seq)
                for rho in sorted(achievable_perms(W, seq)):
                    ninst += 1
                    nr = not routable(W, seq, rho)
                    nnr += nr
                    t = scalar_min_changes(W, seq, rho)
                    if t < W:
                        print('   THEOREM-A VIOLATION', W, seq, rho, t, flush=True); ok = False
                    tmin = min(tmin, t)
                    mb = pattern_bounds(W, n, R, rho, True, 10 ** 7)
                    Mlo, Mhi = min(Mlo, mb[0] / W), min(Mhi, mb[1] / W)
                    xb = pattern_bounds(W, n, R, rho, False, XCAP)
                    if xb is None:
                        Xskip += 1
                    else:
                        Xlo, Xhi = min(Xlo, xb[0] / W), min(Xhi, xb[1] / W)
            print(f'W={W} n={n}: seqs={nseq} (seq,rho)={ninst} nonroutable={nnr} | scalar-any-climb min/W={tmin/W:.3f}'
                  f' | monotone 2-class s/Wm in [{Mlo:.3f},{Mhi:.3f}] | any 2-class s/Wm in [{Xlo:.3f},{Xhi:.3f}]'
                  f' (skipped {Xskip}) [{time.time()-t0:.1f}s]', flush=True)
            ok &= nnr == 0 and min(Mlo, Xlo) >= 1 - 1e-9
    print('PASS_EXPQQ_ENUM' if ok else 'FAIL_EXPQQ_ENUM')


if __name__ == '__main__':
    main()
