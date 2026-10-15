"""Case III (|S| = 1): g_s = (ka_s r_s, -2 al_s) used, e9 not used, no h0-generic line (that case: hz.py, T_Y <= 3).
credit = 1 + dim(D_0 cap q_s^perp) + sum_{k != 0} dim D_k   (D_0 = Def of g_s's partner g'_0).
Every other Y-line pairs 0 with g_s: lies in H* = {t ka_s = al_s gamma}. Pencil lines in H*: only l* (level lam* = -al_s/ka_s,
never vanishing at s since 2 al al' + ka ka' = 1). Generic F_b in H*  <=>  lam* = (3 v + L_b)/4, v = level of M_b^-;
hence such b's with the same M-block share one level (one block), which is neither the M-block nor the l* block, and the
M-block is not the l* block. D_0 cap q_s^perp = D_0 minus 1 dim iff g'_0 is a pencil line whose block contains s.
Labels: 1 = partner pencil block P0 (optional), 2 = l* block (s not in it), 3, 4 = other constant-M blocks, 5 = spare block,
0 = rest.  Usage: s1gs.py [ctrl|noS|sStar]: ctrl drops the H*-level constraints (not binding: 9);
negative controls noS (drop the q_s correction of D_0) and sStar (allow s in the l* block) must reach 10."""
import sys
from ortools.sat.python import cp_model
from s1lib import *
mode = sys.argv[1] if len(sys.argv) > 1 else 'true'
ctrl = mode == 'ctrl'
m = cp_model.CpModel(); K = range(6)
lab = {(c, k): m.NewBoolVar('') for c in C for k in K}
for c in C: m.AddExactlyOne(lab[c, k] for k in K)
if mode != 'sStar': m.Add(lab[s, 2] == 0)
d = {}; used = {}
for k in (1, 2):
    m.Add(sum(lab[c, k] for c in C) <= len(C) - 1)
    inD = {}
    for h, cp in CPART.items():
        v = m.NewBoolVar(''); inD[h] = v
        for c in cp: m.AddImplication(v, lab[c, k])
    e = [m.NewBoolVar('') for j in range(1, 10)]
    for j in range(1, 10):
        m.Add(sum(inD.values()) >= (2 ** j - 1)).OnlyEnforceIf(e[j - 1])
        if j > 1: m.AddImplication(e[j - 1], e[j - 2])
    d[k] = sum(e); used[k] = e[0]
gH = {}; g0 = {}; z = {}
for b in N:
    gH[b] = m.NewBoolVar(''); g0[b] = m.NewBoolVar(''); m.AddBoolOr([gH[b].Not(), g0[b].Not()])
    zs = []
    for k in (1, 2, 3, 4, 5):
        z[b, k] = m.NewBoolVar(''); zs.append(z[b, k])
        for c in Mm[b]: m.AddImplication(z[b, k], lab[c, k])
        m.AddImplication(z[b, k], lab[b, k].Not())
    m.Add(sum(zs) == gH[b] + g0[b])
    if not ctrl:
        m.AddImplication(gH[b], z[b, 2].Not())                 # v != lam*
        m.AddImplication(gH[b], lab[b, 2].Not())               # L_b != lam*
for b in N:
    for b2 in N:
        if b < b2 and not ctrl:
            for k in (1, 2, 3, 4, 5):
                for kk in K:   # same M-block => same level for b, b2
                    m.AddBoolOr([gH[b].Not(), gH[b2].Not(), z[b, k].Not(), z[b2, k].Not(), lab[b, kk].Not(), lab[b2, kk]])
m.Add(used[1] + sum(g0.values()) <= 1)                           # one partner
for b in N:
    if b < star[b]:
        m.Add(gH[b] + gH[star[b]] + g0[b] + g0[star[b]] + used[1] + used[2] <= 3)
m.Add(used[1] + used[2] + sum(gH.values()) + sum(g0.values()) <= 9)
m.Maximize(1 + d[1] - (0 if mode == 'noS' else lab[s, 1]) + d[2] + sum(gH.values()) + sum(g0.values()))
sv = cp_model.CpSolver(); sv.parameters.num_workers = 1; sv.parameters.max_time_in_seconds = 3000
r = sv.Solve(m)
print(mode, sv.StatusName(r), 'max credit', sv.ObjectiveValue(), 'bound', sv.BestObjectiveBound())
if r in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    bl = {k: [c for c in C if sv.Value(lab[c, k])] for k in K}
    print(' block sizes', {k: len(v) for k, v in bl.items()}, 's label', [k for k in K if sv.Value(lab[s, k])],
          'd', sv.Value(d[1]), sv.Value(d[2]), 'gH', [b for b in N if sv.Value(gH[b])], 'g0', [b for b in N if sv.Value(g0[b])])
    print(' recheck dims', dimD(bl[1]), dimD(bl[2]))
