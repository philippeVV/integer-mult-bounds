"""Symmetry-reduced (S_h-equivariant) exact oracle for the bit motif invocation at the upstream scale h = 100.

Same miniature as calib_motif.py (one stage-1 invocation of the 3-subset bit motif, form I - J/9 on Q^h, fixed terminals
X: <t> -> F, Y: 0 -> t^perp, side/centre: 0 -> F), but every gate label is an EQUIVARIANT function of the gate's triple:
  per-triple rows 0,5 (target S) and 2,7 (source T): Stab(S)-invariant subspaces  ell (+) a*V_S (+) b*V_Sc, where
      ell is 0, the trivial 2-plane T2 = <t, j-t>, or one of the lines <t>, <j-t>, <j> and their T2-perps;
      V_S / V_Sc = sum-zero vectors supported on S / S^c (all three summands are mutually orthogonal for I - J/9);
  centre rows 1,3,4,6: S_h-invariant subspaces 0, <j>, j^perp, F.
Edge cost = rank_Q(P_V - P_U) for ANY pair (nested or not).  Minimisation is exact over this family (enumeration of the
4 centre rows x DP over the chain r0 - r2 - r5 - r7).  Ranks: numpy float SVD (tol 1e-7); the recipe's
cost under this model is cross-checked against the exact upstream formula 2v(h-1) + 3ch + vzh.
usage: eqv_motif_y.py h [neg]   (expQT: Y sources may carry any equivariant label l_S, X sinks then l_T^perp)
"""
import sys, itertools
from fractions import Fraction as Fr
from math import comb
import numpy as np

h = int(sys.argv[1]); NEG = len(sys.argv) > 2 and sys.argv[2] == 'neg'
G = np.eye(h) - np.ones((h, h)) / 9.0
S = (0, 1, 2); T = (0, 3, 4)                 # a neighbour pair, |S cap T| = 1
def ind(A): v = np.zeros(h); v[list(A)] = 1; return v
j = np.ones(h)

def sumzero(A):
    A = list(A); out = []
    for k in A[1:]:
        v = np.zeros(h); v[A[0]] = 1; v[k] = -1; out.append(v)
    return out

# exact 2-D trivial part for a triple: basis (t, u = j - t); Gram entries
g_tt = Fr(2); g_tu = Fr(1) - Fr(h, 3); g_uu = Fr(h) - Fr(h * h, 9) - 2 * (3 - Fr(h, 3)) + 2
def ip2(a, b): return a[0] * b[0] * g_tt + (a[0] * b[1] + a[1] * b[0]) * g_tu + a[1] * b[1] * g_uu
def perp2(a):  # (x,y) with ip2((x,y), a) = 0
    c0 = a[0] * g_tt + a[1] * g_tu; c1 = a[0] * g_tu + a[1] * g_uu
    return (c1, -c0) if (c0 or c1) else None
def norm(a):
    if a[0] != 0: return (Fr(1), a[1] / a[0])
    return (Fr(0), Fr(1))
lines = set()
for a in [(Fr(1), Fr(0)), (Fr(0), Fr(1)), (Fr(1), Fr(1))]:
    lines.add(norm(a)); p = perp2(a)
    if p: lines.add(norm(p))
lines = sorted(l for l in lines if ip2(l, l) != 0)     # nondegenerate lines only
triv = ['0', 'T2'] + lines
types = [(tr, a, b) for tr in triv for a in (0, 1) for b in (0, 1)]

def basis_for(tp, A):
    tr, a, b = tp; t = ind(A); u = j - t; B = []
    if tr == 'T2': B += [t, u]
    elif tr != '0': B.append(float(tr[0]) * t + float(tr[1]) * u)
    if a: B += sumzero(A)
    if b: B += sumzero([k for k in range(h) if k not in A])
    return np.array(B) if B else np.zeros((0, h))
central = {'0': np.zeros((0, h)), 'j': j[None, :], 'jp': np.array(sumzero(range(h))), 'F': np.eye(h)}

def proj(B):
    if B.shape[0] == 0: return np.zeros((h, h))
    M = B @ G @ B.T
    return B.T @ np.linalg.solve(M, B @ G)
def rk(M): return int(np.linalg.matrix_rank(M, tol=1e-7))

P_t = {('S', tp): proj(basis_for(tp, S)) for tp in types}
P_t.update({('T', tp): proj(basis_for(tp, T)) for tp in types})
P_c = {k: proj(B) for k, B in central.items()}
P_src_X = proj(ind(T)[None, :])               # X_T source label <t_T>
P_snk_Y = proj(basis_for(((Fr(1), Fr(0)), 0, 0), S)); P_snk_Y = np.eye(h) - P_snk_Y   # t_S^perp
I = np.eye(h); Z = np.zeros((h, h))
d = lambda A, B: rk(A - B)
def dsrc(A):  # from X source
    return rk(A + P_src_X) if NEG else rk(A - P_src_X)

# ---- recipe (upstream labels at stage 1) ----
L_t = ((Fr(1), Fr(0)), 0, 0)
L_tp = None
for tp in types:   # t^perp as a type
    if np.allclose(P_t[('S', tp)], P_snk_Y): L_tp = tp
assert L_tp is not None
L_0, L_F = ('0', 0, 0), ('T2', 1, 1)
rec = dict(r0=L_0, r1='0', r2=L_t, r3='F', r4='0', r5=L_tp, r6='F', r7=L_F)

v = comb(h, 3); z = 3 * comb(h - 3, 2); c = h; W = 2 * v + v * z + c
def total(r):
    X = dsrc(P_t[('T', r['r2'])]) + d(P_t[('T', r['r2'])], P_c[r['r3']]) + d(P_c[r['r3']], P_c[r['r6']]) + \
        d(P_c[r['r6']], P_t[('T', r['r7'])]) + d(P_t[('T', r['r7'])], I)
    Y = d(Z, P_t[('S', r['r0'])]) + d(P_t[('S', r['r0'])], P_c[r['r1']]) + d(P_c[r['r1']], P_c[r['r4']]) + \
        d(P_c[r['r4']], P_t[('S', r['r5'])]) + d(P_t[('S', r['r5'])], P_snk_Y)
    A = d(Z, P_t[('S', r['r0'])]) + d(P_t[('S', r['r0'])], P_t[('T', r['r2'])]) + \
        d(P_t[('T', r['r2'])], P_t[('S', r['r5'])]) + d(P_t[('S', r['r5'])], P_t[('T', r['r7'])]) + d(P_t[('T', r['r7'])], I)
    C = d(Z, P_c[r['r1']]) + d(P_c[r['r1']], P_c[r['r3']]) + d(P_c[r['r3']], P_c[r['r4']]) + \
        d(P_c[r['r4']], P_c[r['r6']]) + d(P_c[r['r6']], I)
    return v * X + v * Y + v * z * A + c * C, (X, Y, A, C)
s_rec_formula = 2 * v * (h - 1) + 3 * c * h + v * z * h + (v if NEG else 0) * 0
s_rec, parts = total(rec)
print(f"h={h} v={v} z={z} W={W} Wm={W*h}  types={len(types)}  recipe via oracle costs = {s_rec} {parts}  "
      f"recipe formula (label sum) = {2*v*(h-1) + 3*c*h + v*z*h}")

# ---- exact minimisation over the family ----
C4 = list(central)
tS = {tp: P_t[('S', tp)] for tp in types}; tT = {tp: P_t[('T', tp)] for tp in types}
# precompute pair tables
dS0 = {a: d(Z, tS[a]) for a in types}
dST = {(a, b): d(tS[a], tT[b]) for a in types for b in types}
dSc = {(a, k): d(tS[a], P_c[k]) for a in types for k in C4}
dTc = {(a, k): d(tT[a], P_c[k]) for a in types for k in C4}
dcc = {(k, l): d(P_c[k], P_c[l]) for k in C4 for l in C4}
dTI = {a: d(tT[a], I) for a in types}
dSY = {a: d(tS[a], P_snk_Y) for a in types}
dXs = {a: dsrc(tT[a]) for a in types}
dc0 = {k: d(Z, P_c[k]) for k in C4}; dcI = {k: d(P_c[k], I) for k in C4}
peryl = {}
YL = [None] + [tp for tp in types if tp != ('0', 0, 0)]
best = None
for yl in YL:
  PyS = tS[yl] if yl else Z; PyT = tT[yl] if yl else Z
  dYs = {a: (rk(tS[a] + PyS) if NEG else rk(tS[a] - PyS)) for a in types}   # Y source label <l_S> (0 = upstream); bit interface carries -P
  dXk = {a: rk(tT[a] - (I - PyT)) for a in types}      # X sink = l'_T^perp (F = upstream)
  for r1, r3, r4, r6 in itertools.product(C4, repeat=4):
      Cc = dc0[r1] + dcc[(r1, r3)] + dcc[(r3, r4)] + dcc[(r4, r6)] + dcI[r6]
      base = c * Cc + v * dcc[(r3, r6)] + v * dcc[(r1, r4)]
      # chain r0 - r2 - r5 - r7 with unary terms
      u0 = {a: v * (dYs[a] + dSc[(a, r1)]) + v * z * dS0[a] for a in types}
      u2 = {a: v * (dXs[a] + dTc[(a, r3)]) for a in types}
      u5 = {a: v * (dSc[(a, r4)] + dSY[a]) for a in types}
      u7 = {a: v * (dTc[(a, r6)] + dXk[a]) + v * z * dTI[a] for a in types}
      # DP
      f2 = {b: min(((u0[a] + v * z * dST[(a, b)], a) for a in types), key=lambda p: p[0]) for b in types}
      f5 = {b: min(((f2[a][0] + u2[a] + v * z * dST[(b, a)], a) for a in types), key=lambda p: p[0]) for b in types}
      f7 = {b: min(((f5[a][0] + u5[a] + v * z * dST[(a, b)], a) for a in types), key=lambda p: p[0]) for b in types}
      val, r7 = min(((f7[b][0] + u7[b], b) for b in types), key=lambda p: p[0])
      r5 = f7[r7][1]; r2 = f5[r5][1]; r0 = f2[r2][1]
      tot = base + val
      peryl[yl] = min(peryl.get(yl, tot), tot)
      if best is None or tot < best[0]:
          best = (tot, dict(yl=yl, r0=r0, r1=r1, r2=r2, r3=r3, r4=r4, r5=r5, r6=r6, r7=r7))
tot, asg = best
for k in YL: print('  Y label', k, 'best s =', peryl[k], ' s - recipe =', peryl[k] - s_rec, ' dim =', (rk(tS[k]) if k else 0))
print(f"equivariant optimum (Y labels allowed) = {tot}  ({'%.6f' % (tot / s_rec)} x recipe)")
print("optimum labels:", asg)
print(f"deficit vs Wm: recipe {W*h - s_rec}, optimum {W*h - tot}")
print(f"RESULT eqv h={h} neg={NEG} recipe={s_rec} opt={tot} Wm={W*h} within5={'yes' if tot <= 1.05*s_rec else 'no'}")
# expQT claim: every nonzero equivariant Y-source label gives s = Wm (loses the whole deficit); the optimum is "no Y label"
ok = asg['yl'] is None and all(peryl[k] == W * h for k in YL if k) and W * h > s_rec
print('PASS_EXPQT_EQV' if ok else 'FAIL_EXPQT_EQV')
