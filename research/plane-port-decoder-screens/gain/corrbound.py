"""EXACT lower bound for 'frame (N_F = G/2) + corrections with denominators dividing 2':
2*N_corr = sum m_k (2-integral terms); mod 2 it equals A_lost (odd G entries on lost pairs) on every required entry,
and B vanishes on free pairs, so by the subadditivity lemma  sum_corr rho_k >= rank_F2(B o (A_lost (x) J)).
Usage: corrbound.py file.json [key]   (planes r from framecorr_*.json / gamma_runs.json / mf files in expG coords)"""
import sys, json, os
from base import *
HERE = os.path.dirname(os.path.abspath(__file__))
emap = {int(k): v for k, v in json.load(open(os.path.join(HERE, '../lower/e8map.json'))).items()}
ip = lambda a, b: sum(x * y for x, y in zip(a, b)) // 4
R = [emap[q] for q in P]; frame = [emap[x] for x in sorted(set(emap) - set(P))]
def bound(r):
    F = free_matrix(r); E = planes(r); rows = []
    lost = 0
    for a in range(n):
        for x in E[a]:
            v = 0; bit = 0
            for b in range(n):
                on = a != b and not F[a][b] and ip(R[a], R[b]) % 2 != 0
                for y in E[b]:
                    if on and dot(x, y): v |= 1 << bit
                    bit += 1
            rows.append(v)
        lost += sum(1 for b in range(n) if a != b and not F[a][b] and ip(R[a], R[b]) % 2)
    D = [[a for a in range(n) if ip(R[a], f)] for f in frame]
    fc = [term_cost(r, d, d) for d in D]
    return lost // 2, sum(fc), rank2(rows)
d = json.load(open(sys.argv[1])); runs = d if isinstance(d, list) else [d]
for run in runs:
    r = run.get('r') or run.get('R')
    lost, fc, rb = bound(r)
    print(sys.argv[1], 'lost pairs', lost, 'frame cost', fc, 'corr bound', rb, '=> frame+half-int corr >=', fc + rb, flush=True)
