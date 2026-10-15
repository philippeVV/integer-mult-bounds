"""|S| = 1 structure library (s = line 0). Y-row model after torus/GL9 gauge (README sec 1):
  a in E = orth(s): row (r_a, 0)          b in N = nonorth(s): r_b = r_s/2 + p_b (sign fixed r_b.r_s = 1), level L_b
  Y-line F = (u0 in r_s^perp, t, gamma): zero at a in E iff r_a.u0 = 0; at b in N iff p_b.u0 + t + gamma L_b = 0;
  at s iff 2t + gamma L_s = 0  (L_s = kap'_s/alp'_s, infinity if alp'_s = 0)."""
import sys, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from e8 import *
from maxdef import span2
R = np.array(ROOT, float) / 2
s = 0
E = [a for a in range(NL) if IP[a][s] == 0 and a != s]; N = [a for a in range(NL) if IP[a][s] != 0 and a != s]
C = [s] + N
sg = {b: np.sign(R[b] @ R[s]) for b in N}
p = {b: sg[b] * R[b] - R[s] / 2 for b in N}
PP = {(b, c): round(float(p[b] @ p[c]) * 2) for b in N for c in N}     # 2 p_b.p_c in {3,1,-1,-3}
star = {b: next(c for c in N if PP[b, c] == -3) for b in N}
Mm = {b: [c for c in N if PP[b, c] == -1] for b in N}                  # M_b^-  (27)
Mp = {b: [c for c in N if PP[b, c] == 1] for b in N}                   # M_b^+  (27)
CL = {h: [a for a in range(NL) if dotq(LAB[a], h)] for h in range(1, 511)}
CPART = {h: frozenset(a for a in O if a in C) for h, O in CL.items()}  # O_h cap C(s)
H0 = next(h for h, O in CL.items() if sorted(O) == sorted(N))
HB = {}
for h, O in CL.items():
    OE = [a for a in O if a in E]
    if len(OE) == 36:
        b = [b for b in N if b in O and set(Mm[b]) <= set(O)]
        assert len(b) == 1; HB[b[0]] = h
assert len(HB) == 56
for b in N: assert HB[b] ^ HB[star[b]] == H0
def dimD(P):
    """pencil Def dim for a level block P (subset of C)."""
    P = set(P); return len([h for h in CL if CPART[h] <= P]).bit_length()
