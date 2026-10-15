"""expBF base: ports from pinned bundle (read-only import of lower/ports.py), instances."""
import sys, os, json
RES = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RES, 'lower'))
from ports import P, dot, rchoices
n = len(P)
def load_r(name):
    f = {'a23': 'construct2/a23_first.json', 'a24': 'construct2/a24_first.json', 'amin25': 'construct2/amin25.json',
         'amin_s1': 'construct2/amin_s1.json', 'rW': 'gain/rW_1.json'}[name]
    d = json.load(open(os.path.join(RES, f)))
    if isinstance(d, list): d = d[0]
    return d['r']
def planes(r): return [(P[i], r[i]) for i in range(n)]
def orth(Ea, Eb): return all(dot(x, y) == 0 for x in Ea for y in Eb)
def free(r):
    E = planes(r); return [[i != j and orth(E[i], E[j]) for j in range(n)] for i in range(n)]
def reduce(vs):
    basis = []
    for v in vs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return basis
def rank2(rows): return len(reduce(rows))
def pair_rank(U, V):
    bu, bv = reduce(U), reduce(V)
    return rank2([sum(dot(a, b) << j for j, b in enumerate(bv)) for a in bu])
def cost(r, R, A):
    return pair_rank([w for i in R for w in (P[i], r[i])], [w for j in A for w in (P[j], r[j])])
