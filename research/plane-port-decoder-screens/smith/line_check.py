"""Smith lemma in the LINE model on the incumbent G/2 (8 frame terms (1/4) c_k c_k^T): 112 <= sum_k r_{d_k}(Bt_line[D_k,D_k])."""
import json, os
from fractions import Fraction as Fr
from common import *
from smithnp import r_j
from smith import v2, g
HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gain')
emap = {int(k): v for k, v in json.load(open(os.path.join(HERE, '../lower/e8map.json'))).items()}
ip = lambda a, b: sum(x * y for x, y in zip(a, b)) // 4
Rt = [emap[q] for q in P]; frame = [emap[x] for x in sorted(set(emap) - set(P))]
c = [[ip(Rt[a], f) for f in frame] for a in range(n)]
terms = [({a: Fr(c[a][k], 4) for a in range(n) if c[a][k]}, {a: Fr(c[a][k]) for a in range(n) if c[a][k]}) for k in range(8)]
N = [[sum(d.get(a, 0) * t.get(b, 0) for d, t in terms) for b in range(n)] for a in range(n)]
assert all(N[a][b] == (1 if a == b else 0) for a in range(n) for b in range(n) if a == b or dot(P[a], P[b]) == 1)
tot = 0
for d, t in terms:
    e = min(v2(x) for x in d.values()) + min(v2(y) for y in t.values()); R = sorted(d); A = sorted(t)
    M = [[dot(P[a], P[b]) for b in A] for a in R]
    rho = r_j(M, 0); h = r_j(M, -e); tot += h
    print('term |R|', len(R), 'rho', rho, 'depth', -e, 'r_depth', h, 'g', g(-e, rho))
print('LINE incumbent exact OK; C =', sum(r_j([[dot(P[a], P[b]) for b in sorted(t)] for a in sorted(d)], 0) for d, t in terms), ' Smith sum', tot, '>= 112:', tot >= 112)
