"""EXACT. Gain inequality: split terms into K+ (term matrix 2-integral) and K- (rest). Then
   C >= NAIVE - gain(K-),  gain(K-) = rank_F2( B o (Nbar_- (x) J) ) - sum_{K-} rho_k,
where N_- = sum_{K-} terms is 2-integral on required entries (= I - N_+ there), Nbar_- its reduction on required entries,
NAIVE = rank of blockdiag (224 planes, 112 lines).  [lemma on K+, then rank(D - Y) >= rank D - rank Y]
Checks: (1) LINE incumbent (frame decoder N = (1/4) sum c c^T): exact identity, gain = 112 - 64 = 48 => bound 64 tight.
        (2) plane decoders from rectF_*.json (all terms non-integral): gain = 224 - C."""
import sys, json, os
from fractions import Fraction as Fr
from base import *
HERE = os.path.dirname(os.path.abspath(__file__))
emap = {int(k): v for k, v in json.load(open(os.path.join(HERE, '../lower/e8map.json'))).items()}
ip = lambda a, b: sum(x * y for x, y in zip(a, b)) // 4
Rt = [emap[q] for q in P]; frame = [emap[x] for x in sorted(set(emap) - set(P))]
c = [[ip(Rt[a], f) for f in frame] for a in range(n)]
D = [[a for a in range(n) if c[a][k]] for k in range(8)]
def v2int(x): return Fr(x).denominator % 2 == 1
def gain(terms, basis, required):
    """basis[a] = list of F2 vectors (label basis); required[a][b] bool; terms (dict,dict) Fractions."""
    Km = [(d, t) for d, t in terms if not all(v2int(Fr(x) * Fr(y)) for x in d.values() for y in t.values())]
    Nm = {}
    for d, t in Km:
        for i, x in d.items():
            for j, y in t.items(): Nm[(i, j)] = Nm.get((i, j), 0) + Fr(x) * Fr(y)
    rows = []
    for a in range(n):
        for x in basis[a]:
            v = 0; bit = 0
            for b in range(n):
                val = Nm.get((a, b), 0)
                on = False
                if required[a][b]:
                    assert v2int(val), 'N_- not 2-integral on a required entry'
                    on = (val.numerator % 2 == 1)
                for y in basis[b]:
                    if on and dot(x, y): v |= 1 << bit
                    bit += 1
            rows.append(v)
    cost = sum(pair_rank([w for i in d for w in basis[i]], [w for j in t for w in basis[j]]) for d, t in Km)
    return len(Km), rank2(rows), cost
# (1) line incumbent
lt = [({a: Fr(c[a][k], 4) for a in D[k]}, {a: Fr(c[a][k]) for a in D[k]}) for k in range(8)]
req = [[a == b or dot(P[a], P[b]) == 1 for b in range(n)] for a in range(n)]
N = [[sum(Fr(d.get(a, 0)) * Fr(t.get(b, 0)) for d, t in lt) for b in range(n)] for a in range(n)]
assert all(N[a][b] == (1 if a == b else 0) for a in range(n) for b in range(n) if req[a][b])
k, rk, cost = gain(lt, [[q] for q in P], req)
print('LINE incumbent: exact OK; non-integral terms', k, 'rank Y', rk, 'cost', cost, 'gain', rk - cost, '=> C >=', 112 - (rk - cost))
for f in sys.argv[1:]:
    d = json.load(open(f)); r = d['r']
    terms = [({int(i): Fr(x) for i, x in dd.items()}, {int(j): Fr(y) for j, y in tt.items()}) for dd, tt in d['terms']]
    F = free_matrix(r); reqp = [[not F[a][b] for b in range(n)] for a in range(n)]
    k, rk, cost = gain(terms, [[P[a], r[a]] for a in range(n)], reqp)
    print(f, 'C', d['C'], 'non-integral terms', k, 'rank Y', rk, 'cost', cost, 'gain', rk - cost)
