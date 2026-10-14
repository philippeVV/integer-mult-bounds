"""max #8-cliques touched by a rectangle (R,A) with linerank<=rho (rho=1: also R∩A=∅)."""
import sys
from common import *
from ortools.sat.python import cp_model
rho = int(sys.argv[1]); tl = float(sys.argv[2]) if len(sys.argv) > 2 else 300
m = cp_model.CpModel()
x = [m.NewBoolVar(f'x{a}') for a in range(n)]; y = [m.NewBoolVar(f'y{a}') for a in range(n)]
f = [[m.NewBoolVar('') for _ in range(9)] for _ in range(rho)]; g = [[m.NewBoolVar('') for _ in range(9)] for _ in range(rho)]
def lin(vec, q):
    v = m.NewBoolVar(''); bits = [vec[i] for i in range(9) if q >> i & 1]
    m.AddBoolXOr(bits + [v.Not()]) if bits else m.Add(v == 0)  # v = xor(bits)
    return v
s = [[lin(f[l], P[a]) for l in range(rho)] for a in range(n)]
t = [[lin(g[l], P[b]) for l in range(rho)] for b in range(n)]
for a in range(n):
    if rho == 1: m.AddBoolOr([x[a].Not(), y[a].Not()])
    for b in range(n):
        prods = []
        for l in range(rho):
            p = m.NewBoolVar(''); m.AddMultiplicationEquality(p, [s[a][l], t[b][l]]); prods.append(p)
        # xor(prods) == q_a.q_b  when x_a & y_b
        if dot(P[a], P[b]): m.AddBoolXOr(prods + [m.NewConstant(1)]).OnlyEnforceIf([x[a], y[b]]) if False else m.AddBoolXOr(prods).OnlyEnforceIf([x[a], y[b]])
        else: m.AddBoolXOr(prods + [m.NewConstant(1)]).OnlyEnforceIf([x[a], y[b]])
z = []
for S in CL8:
    zz = m.NewBoolVar(''); m.AddBoolOr([x[a] for a in S]).OnlyEnforceIf(zz); m.AddBoolOr([y[a] for a in S]).OnlyEnforceIf(zz); z.append(zz)
m.Maximize(sum(z))
sv = cp_model.CpSolver(); sv.parameters.max_time_in_seconds = tl; sv.parameters.num_workers = 2
st = sv.Solve(m)
R = [a for a in range(n) if sv.Value(x[a])]; A = [a for a in range(n) if sv.Value(y[a])]
tc = sum(1 for S in CL8 if set(S) & set(R) and set(S) & set(A))
print('rho', rho, sv.StatusName(st), 'best', sv.ObjectiveValue(), 'bound', sv.BestObjectiveBound(), 'recheck touch', tc, 'linerank', linerank(R, A), '|R|', len(R), '|A|', len(A), 'RA', len(set(R)&set(A)), 'of', len(CL8), flush=True)
