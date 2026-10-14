"""Lemma C' test: cross-class constraints + 'each class is served by >= k rank-1 terms'
(term x serves class a iff u_x and v_x are both not identically zero on class a)."""
import sys
from pycryptosat import Solver
from sat1 import build
from cases import *
k = int(sys.argv[1])
for ci in range(len(OPEN)):
    cl_, n, nc, nt, var = build(ci)
    A, c = OPEN[ci]; two, one = split(A, c); E6, cl = classes(two)
    S = Solver(); S.add_clauses(cl_); top = max(abs(l) for cc in cl_ for l in cc) + 10000
    srv = {}
    for x in one:
        for a in two:
            W = [r for r in cl[a] if r not in F0 and GM[r][x] == 0]
            top += 1; s = top; srv[x, a] = s
            S.add_clause([-s] + [var['U', x, r] for r in W]); S.add_clause([-s] + [var['V', x, r] for r in W])
    for a in two:
        xs = [srv[x, a] for x in one]
        if k == 2:
            for i in range(3):   # at least 2 of 3: every pair of negations fails -> each 2-subset of complement
                S.add_clause([xs[j] for j in range(3) if j != i])
        elif k == 1:
            S.add_clause(xs)
    sat, sol = S.solve()
    print(ci, OPEN[ci], 'k', k, 'SAT' if sat else 'UNSAT', flush=True)
