"""expNM exact decoders for the field-line families (exact integer/Fraction arithmetic, F2 pairing ranks).

D0  identity decoder: V terms, each costs u -> gamma = 1.
D1  free-clique (frame) decoder: partition the ports into pairwise-free sets ('self-polar frames') found by CP-SAT,
    N = sum_C 1_C 1_C^T (+ identity on leftovers).  N_aa = 1, off-diagonal nonzeros only inside free sets -> valid.
    It realises the LM minrank V/omega; each frame term spans all of F_q^n so costs h -> gamma = 1 exactly
    (the u-dim analogue of PR2244 Theorem C).
usage: uv run -q --no-project --python 3.12 --with numpy --with ortools python -I -B exact.py u n dot|herm [plus1]"""
import os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fam import build, f2_label, f2_form
from ortools.sat.python import cp_model

u, n, kind = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
plus1 = 'plus1' in sys.argv
F, B, P = build(u, n, kind, plus1)
V = len(P); h = u * n
G = f2_form(F, n, kind)
Gm = [sum(G[i][j] << j for j in range(h)) for i in range(h)]   # row bitmasks of the F2 Gram


def ip(x, y):   # F2 form on bitvectors
    t = 0
    for i in range(h):
        if x >> i & 1: t ^= bin(Gm[i] & y).count('1') & 1
    return t


def basis(vs):
    b = []
    for r in vs:
        for x in b: r = min(r, r ^ x)
        if r: b.append(r)
    return b


def rank2(rows):
    return len(basis(rows))


def pair_rank(R, A):
    Rb, Ab = basis(R), basis(A)
    return rank2([sum(ip(x, y) << j for j, y in enumerate(Ab)) for x in Rb])


L = [f2_label(F, x) for x in P]
for a in range(V): assert pair_rank(L[a], L[a]) == u, 'label degenerate'
free = np.array([[a != b and B(P[a], P[b]) == 0 for b in range(V)] for a in range(V)])
# cross-check: F2 pairing rank is u exactly on non-free pairs and 0 on free pairs (field-line alphabet)
import random
rnd = random.Random(1)
for _ in range(300):
    a, b = rnd.randrange(V), rnd.randrange(V)
    if a != b: assert pair_rank(L[a], L[b]) == (0 if free[a, b] else u)
print(f"family u={u} n={n} {kind}{' +1' if plus1 else ''}: V={V}, h={h}, field-line PR in {{0,{u}}} (checked on 300 pairs)")

# D0
print(f"D0 identity: C = {u * V}, gamma = 1")

# D1 maximum packing of pairwise-free sets of size omega = n (CP-SAT over all free n-cliques)
import itertools
adj = [set(np.nonzero(free[a])[0].tolist()) for a in range(V)]
cl = []
def ext(c, cand):
    if len(c) == n: cl.append(tuple(c)); return
    for b in sorted(cand):
        if b > c[-1]: ext(c + [b], cand & adj[b])
for a in range(V): ext([a], adj[a])
m = cp_model.CpModel(); y = [m.NewBoolVar('') for _ in cl]
for a in range(V): m.Add(sum(y[i] for i, c in enumerate(cl) if a in c) <= 1)
m.Maximize(sum(y)); s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = 120; s.parameters.num_workers = 2
s.Solve(m)
chosen = [cl[i] for i in range(len(cl)) if s.Value(y[i])]
cover = {a for c in chosen for a in c}; left = [a for a in range(V) if a not in cover]
terms = [(list(c), list(c)) for c in chosen] + [([a], [a]) for a in left]
N = {}
for R, A in terms:
    for a in R:
        for b in A: N[a, b] = N.get((a, b), 0) + 1
ok = all(N.get((a, a), 0) == 1 for a in range(V)) and all(free[a, b] for (a, b), v in N.items() if a != b and v)
costs = [pair_rank([x for a in R for x in L[a]], [x for b in A for x in L[b]]) for R, A in terms]
C = sum(costs)
print(f"D1 frame decoder: {len(cl)} free {n}-cliques, packing {len(chosen)} frames + {len(left)} singletons = {len(terms)} terms "
      f"(LM minrank ~ V/omega = {V / n:.1f}); exact valid: {ok}; term costs {sorted(set(costs))}; C = {C}, gamma = {Fr(C, u * V)}")
assert ok
print('PASS_EXPNM_EXACT')
