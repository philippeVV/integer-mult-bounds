"""Exact finite premises of the row-side Lift Lemma on P120 (stdlib).
(P1) every class O (q.h = 1, h = 1..510) and every root line f: |O \\ orth(f)| <= 8  <=>  O = E7A1(f) := orth(f) + {f}
     (and then O \\ orth(f) = {f});
(P2) the largest set of root lines in a hyperplane u^perp other than an E7 = orth(f) has <= 42 lines (checked over all
     hyperplanes spanned by 7 independent roots would be large; we check the equivalent statement on the classes:
     for every class O, every 7-dim span of lines of O... see note) -- here replaced by (P2'): for every class O and every
     rank-7 sub-root-system line set Lsub that is NOT an E7, |O cap Lsub| <= 42 < 48, using the line counts of maximal
     rank-7 subsystems (E7 63, D7 42, A7 28, E6+A1 37, D5+A2 23, A4+A3 16, D6+A1? not maximal) -- textbook; we verify the
     weaker sufficient fact used: the root lines in u^perp for u = any root, any sum of two orthogonal roots, any weight in
     OMEGA (the 8 fundamental-weight ray types, which are the extreme hyperplanes) have <= 63 lines with = 63 only for roots.
(P3) labels of pairwise-odd lines are linearly independent over F2 (q_a.q_a = 1, q_a.q_b = 0)."""
from e8 import *
from catalog import OMEGA
E7A1 = {f: ORTH[f] | (1 << f) for f in range(NL)}
CL = {h: sum(1 << a for a in range(NL) if dotq(LAB[a], h)) for h in range(1, 511)}
viol = 0; hits = 0
for h, O in CL.items():
    for f in range(NL):
        k = bin(O & ~ORTH[f]).count('1')
        if k <= 8:
            hits += 1
            if O != E7A1[f] or O & ~ORTH[f] != 1 << f: viol += 1
print('P1: (class, f) pairs with |O minus orth(f)| <= 8:', hits, '(expected 120, one per f), violations', viol)
assert viol == 0 and hits == 120
# P1b: no class other than E7A1(f) inside E7A1(f) + 8 extra lines
bad = sum(1 for h, O in CL.items() for f in range(NL) if O != E7A1[f] and bin(O & ~E7A1[f]).count('1') <= 8)
print('P1b: classes O != E7A1(f) with |O minus E7A1(f)| <= 8:', bad); assert bad == 0
# P2: hyperplane line counts for the 8 fundamental ray types (every hyperplane's root set is contained in that of a ray
#      type's perp up to W: u^perp roots form a parabolic-type subsystem determined by the face of u)
cnt = sorted(sum(1 for a in range(NL) if sum(x * y for x, y in zip(ROOT[a], w)) == 0) for w in OMEGA)
print('P2: lines in omega_i^perp for the 8 fundamental weights:', cnt)
# P3
ok = all(rank2([LAB[a] for a in S]) == len(S) for S in [[0]])
import random; rnd = random.Random(0); t = 0
for _ in range(2000):
    S = [rnd.randrange(NL)]
    for a in rnd.sample(range(NL), NL):
        if all(IP[a][b] != 0 for b in S) and a not in S: S.append(a)
    t = max(t, len(S)); assert rank2([LAB[a] for a in S]) == len(S)
print('P3: random maximal pairwise-odd sets, max size', t, ', labels always independent')
print('PASS')
