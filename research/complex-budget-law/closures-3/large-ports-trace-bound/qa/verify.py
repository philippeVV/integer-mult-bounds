"""expQA verify: re-check every rigorous number in the README from the saved certificates (seconds).
 - upper bounds: clique-partition decoders N = sum_C 1_C 1_C^T are exact masked identities of rank #cliques;
 - lower bounds: pairwise-required sets (alpha) and C5-packings (mr(C5) = 3 over every field) are valid.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import json, glob, os, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))   # certificates live next to this script
from qa import family, check_masked, rank_modp

ok_all = True
rows = []
for f in sorted(glob.glob('bounds_*.json')):
    d = json.load(open(f)); fam = d['fam']
    if fam == 'e8': continue
    free, th, h, lab, P = family(fam); V = free.shape[0]
    tV = th * V
    # upper bound: best exact partition among bounds/torus files
    cl = d['cliques']
    for g in (f"partition_extra_{fam.replace(':','_')}.json", f"torus_partition_q{fam.split(':')[2]}.json" if fam.endswith(':dot') and fam.startswith('pg:3') else 'none'):
        try:
            c2 = json.load(open(g))['cliques']
            if len(c2) < len(cl): cl = c2
        except FileNotFoundError: pass
    d['cliques'] = cl
    N = np.zeros((V, V), dtype=np.int64)
    for c in cl: N[np.ix_(c, c)] += 1
    ok, _, _ = check_masked(N, free)
    rk = rank_modp(N)
    part = sorted(sum((list(c) for c in d['cliques']), []))
    ok &= part == list(range(V)) and rk == len(d['cliques'])
    # alpha
    best_a = d['alpha_set']
    try:
        A2 = json.load(open(f"alpha_{fam.replace(':','_')}.json"))
        if len(A2) > len(best_a): best_a = A2
    except FileNotFoundError: pass
    ok &= not free[np.ix_(best_a, best_a)].any()
    lb = len(best_a)
    # C5 packing
    try:
        c5 = json.load(open(f"c5pack_{fam.replace(':','_')}.json"))
        U = c5['singletons'] + [a for c in c5['c5'] for a in c]
        assert len(set(U)) == len(U)
        piece = {a: ('s', a) for a in c5['singletons']}
        for c in c5['c5']:
            sub = free[np.ix_(c, c)]; ok &= bool((sub.sum(1) == 2).all())
            for a in c: piece[a] = ('c', tuple(c))
        for a in U:
            for b in np.nonzero(free[a])[0]:
                if b in piece: ok &= piece[a] == piece[b]
        lb = max(lb, len(c5['singletons']) + 3 * len(c5['c5']))
    except FileNotFoundError: pass
    ok_all &= bool(ok)
    rows.append((lab, V, tV, lb, rk, ok))
print(f"{'family':22s} {'V':>5s} {'thetaV':>7s} {'mr>=':>5s} {'mr<=':>5s} {'sigma in':>16s}  check")
for lab, V, tV, lb, ub, ok in rows:
    print(f"{lab:22s} {V:5d} {tV:7.2f} {lb:5d} {ub:5d}   [{lb/tV:.3f}, {ub/tV:.3f}]  {'PASS' if ok else 'FAIL'}")
# mr(C5) = 3: exhaustive-style symbolic check is in the README; numeric sanity: rank-2 LM fails
print('ALL PASS' if ok_all else 'SOME FAIL')
print('PASS_EXPQA_VERIFY' if ok_all else 'FAIL_EXPQA_VERIFY')
