"""expJ base: ports (expG model: std dot on F2^9), plane spans, F2 term cost, exact decoder verification."""
import sys, os, json
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lower'))
from ports import P, dot, rchoices
n = len(P)
NB = [[j for j in range(n) if j != i and dot(P[i], P[j]) == 0] for i in range(n)]
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
def form_rank(U): return pair_rank(U, U)
def planes(r): return [(P[i], r[i]) for i in range(n)]
def orth_pl(Ea, Eb): return all(dot(x, y) == 0 for x in Ea for y in Eb)
def free_matrix(r):
    E = planes(r)
    return [[i != j and orth_pl(E[i], E[j]) for j in range(n)] for i in range(n)]
def term_cost(r, rows, cols):
    return pair_rank([w for i in rows for w in (P[i], r[i])], [w for j in cols for w in (P[j], r[j])])
def verify(r, terms):
    """terms: list of (dict row->Fraction, dict col->Fraction). Exact: N_aa=1, N_ab=0 on non-orth pairs.
    returns (bad_count, C, costs)."""
    assert all(dot(r[i], P[i]) == 0 and dot(r[i], r[i]) == 1 for i in range(n))
    F = free_matrix(r)
    N = [[Fr(0)] * n for _ in range(n)]
    for d, t in terms:
        for i, x in d.items():
            for j, y in t.items(): N[i][j] += Fr(x) * Fr(y)
    bad = sum(1 for i in range(n) for j in range(n) if not F[i][j] and N[i][j] != (1 if i == j else 0))
    costs = [term_cost(r, list(d), list(t)) for d, t in terms]
    return bad, sum(costs), costs
