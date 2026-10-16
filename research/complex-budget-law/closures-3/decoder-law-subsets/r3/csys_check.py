"""c-system dimension (extra = dim - 8 - |S|; 0 = concurrent only) on structured and random S, with controls."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys, random, itertools
from r3lib import *
rnd = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
def span_lines(S):
    B = RF[list(S)]; u, s, vt = np.linalg.svd(B); Q = vt[(s > 1e-9).sum() * 0:(s > 1e-9).sum()]
    res = (RF * RF).sum(1) - ((RF @ Q.T) ** 2).sum(1); return [a for a in range(NL) if res[a] < 1e-9]
def report(name, S):
    S = sorted(S); bl = bad_lines(S)
    e_all = csys(S); e_T = csys(S, [a for a in range(NL) if a not in S and a not in bl])
    print(f'{name:34s} |S|={len(S):3d} rank={rank_S(S)} clique={is_clique(S)!s:5s} max|S∩orth|={max(sum(1 for s in S if IP[a][s]==0 and s!=a) for a in range(NL)):3d} bad={len(bl):3d} extra(all good lines)={e_all} extra(T-lines only)={e_T}', flush=True)
    return e_all
a3 = span_lines([0, 1, 2]) if len(span_lines([0, 1, 2])) == 6 else None
# find an A3 flat
for tri in itertools.combinations(range(30), 3):
    L = span_lines(tri)
    if len(L) == 6: A3 = L; break
for tri in itertools.combinations(range(40), 3):
    L = span_lines(tri)
    if len(L) == 4: A2A1 = L; break
fl4 = None
for q in itertools.combinations(range(40), 4):
    L = span_lines(q)
    if len(L) == 12 and rank_S(L) == 4: D4 = L; break
print('A3 lines', A3, 'A2+A1', A2A1, 'D4', D4)
A2 = [s for s in A2A1 if sum(1 for t in A2A1 if IP[s][t] == 0) == 1]
ok = True
ok &= report('A2 triple (control: expect 1)', A2) == 1
ok &= report('A2 + orth line', A2A1) == 0
for q in itertools.combinations(A3, 4):
    if rank_S(q) == 3: ok &= report('4 lines of A3 ' + ('clique' if is_clique(q) else 'non-clique'), q) == 0
ok &= report('A3 (6 lines)', A3) == 0
ok &= report('D4 (12 lines)', D4) == 0
for k in [5, 6, 8, 9, 10, 12]:
    for q in itertools.combinations(D4, k):
        if not is_clique(q): ok &= report(f'{k} lines of D4', q) == 0; break
# random dependent non-clique S of sizes 4..26, biased: random S inside a random rank-r flat + extras
for trial in range(60):
    k = rnd.randint(4, 26)
    base = span_lines(rnd.sample(range(NL), rnd.randint(3, 6)))
    S = set(rnd.sample(base, min(len(base), rnd.randint(3, k))))
    while len(S) < k: S.add(rnd.randrange(NL))
    if is_clique(S) or rank_S(S) == len(S): continue
    ok &= report(f'random #{trial}', S) == 0
# large S with bad lines: the E7-cocircuit Q = orth(a) minus E6 (27 lines) and E6 itself (36)
a = 0; E7 = ORTHL[a]; s1 = next(b for b in range(NL) if IP[a][b] != 0 and b != a)
E6 = [b for b in E7 if IP[s1][b] == 0]; Q = [b for b in E7 if b not in E6]
report('Q = orth(a) minus E6 (27, R3-residual)', Q)
report('E6 lines (36)', E6)
report('control: A2 triple, T-only, no lines dropped', A2)
# negative control: drop the constraints of a D6 cocircuit (10 lines) for an orthogonal pair in S -> extra freedom appears
S = [b for b in A3]; s, t = next((x, y) for x in S for y in S if x < y and IP[x][y] == 0)
D6 = [b for b in range(NL) if IP[s][b] == 0 and IP[t][b] == 0 and b not in (s, t)]
# a 20-line hyperplane flat of D6
import itertools as it
for q in it.combinations(D6, 5):
    if rank_S(q) == 5:
        H = [b for b in D6 if rank_S(list(q) + [b]) == 5]
        if len(H) == 20: break
drop = set(D6) - set(H)
e0 = csys([s, t], [x for x in range(NL) if x not in drop]); e = csys(S, [x for x in range(NL) if x not in drop])
print('NEG CONTROL S={s,t} orthogonal pair, the 10 D6-cocircuit lines dropped: extra =', e0, 'PASS' if e0 > 0 else 'FAIL', '| same drop inside A3 (other pairs re-force concurrency): extra =', e)
e = e0
ok &= e > 0
print('ALL concurrent as predicted' if ok else 'MISMATCH')
print('PASS_EXPR3_CSYS' if ok else 'FAIL_EXPR3_CSYS')
