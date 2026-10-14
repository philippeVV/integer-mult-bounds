# Exact replay of the selected canonical word: split y += J x into drop-carried (via cut) and lossless parts.
import sys, json, pickle
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lower'))
from fractions import Fraction as Q
from ports import P, W, dot
v, R = W['v'], W['R']; n = 2*v + R
def apply(e, rows, cutmode=False):
    if e[0] == 'scale':
        _, t, a, b, *_ = e; q = Q(a, b)
        rows[t] = {j: q*x for j, x in rows[t].items()}; return
    _, t, s, a, b, *_ = e; q = Q(a, b)
    src = {('c', s): Q(1)} if cutmode else rows[s]
    for j, x in list(src.items()):
        z = rows[t].get(j, Q(0)) + q*x
        if z: rows[t][j] = z
        else: rows[t].pop(j, None)
rows = [{r: Q(1)} for r in range(n)]
# compensation as in hpe_audit_and_transport: built from dirty part of main
main = W['A'] + W['cut'] + W['B']
for e in W['precondition'] + W['compensation'] + W['A']: apply(e, rows)
cent = W['centers']
assert all(e[2] in cent for e in W['cut']), 'cut src not centre'
assert all(v <= e[1] < 2*v for e in W['cut']), 'cut dst not target'
D = {c: {j: x for j, x in rows[c].items() if isinstance(j, int) and j < v} for c in cent}
Ddirty = {c: {j: x for j, x in rows[c].items() if not (isinstance(j, int) and j < v)} for c in cent}
for e in W['cut']: apply(e, rows, cutmode=True)
for e in W['B'] + W['cleanup'] + W['tail']: apply(e, rows)
# final y rows
ok = True
L = [[Q(0)]*v for _ in range(v)]; T = [[Q(0)]*8 for _ in range(v)]
for t in range(v):
    row = rows[v+t]
    for j, x in row.items():
        if isinstance(j, tuple): T[t][cent.index(j[1])] = x
        elif j < v: L[t][j] = x
        elif j == v+t: assert x == 1
        else: ok = False
print('other deps in y rows?', not ok)
Dm = [[D[c].get(a, Q(0)) for a in range(v)] for c in cent]
N = [[sum(T[t][k]*Dm[k][a] for k in range(8)) for a in range(v)] for t in range(v)]
J = [[N[t][a] + L[t][a] for a in range(v)] for t in range(v)]
print('J == I:', all(J[t][a] == (1 if t == a else 0) for t in range(v) for a in range(v)))
# classify by label parity
from collections import Counter
cN = Counter(); cL = Counter()
for t in range(v):
    for a in range(v):
        kind = 'diag' if t == a else ('even(q.q=1)' if dot(P[t], P[a]) else 'odd(q.q=0)')
        cN[(kind, N[t][a])] += 1; cL[(kind, L[t][a])] += 1
print('N values by class:', sorted(cN.items(), key=str))
print('L values by class:', sorted(cL.items(), key=str))
print('dirty content of centres at cut (nonzero?):', {c: len(Ddirty[c]) for c in cent})
nP = len(P)
ok = all((1 + 2*N[a][b]) % 2 == dot(P[a], P[b]) for a in range(nP) for b in range(nP))
print('label Gram B0 = (J + 2N) mod 2 entrywise:', ok, '| N integral (B0 a possible target)?', all(N[a][b].denominator == 1 for a in range(nP) for b in range(nP)))
if len(sys.argv) > 1: pickle.dump(dict(N=N, L=L, T=T, D=Dm), open(sys.argv[1], 'wb'))
