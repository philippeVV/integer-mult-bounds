"""Exact min of the phase-frame cost s on Braverman-He scenario A (CP-SAT, 2 workers).
usage: solve.py m {cost|rank} timelimit_s [extra_constraint]
Terminal frames are eliminated exactly (min over beta_in of c(F - beta_in) + c(beta_in + I - L) = c(F + L + I)),
gauge beta[0] = 0.  Wm = 157 m; flow bound s >= 147 m."""
import sys, time
from ortools.sat.python import cp_model
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
import beta as B

m = int(sys.argv[1]); var = sys.argv[2]; TL = float(sys.argv[3])
nb, rk, cost, I = B.tables(m)
tab = cost if var == 'cost' else rk
md = cp_model.CpModel()
bits = [[md.NewBoolVar(f'b{v}_{k}') for k in range(nb)] for v in range(B.NV)]
for k in range(nb): md.Add(bits[0][k] == 0)
terms = []
allowed = [tuple(((x >> k) & 1) for k in range(nb)) + (tab[x],) for x in range(1 << nb)]
def edge(u, v, off):
    d = []
    for k in range(nb):
        dk = md.NewBoolVar('')
        # dk = bits[u][k] xor bits[v][k] xor off_k
        md.AddBoolXOr([bits[u][k], bits[v][k], dk.Not() if not (off >> k & 1) else dk])
        d.append(dk)
    c = md.NewIntVar(0, m + 1, '')
    md.AddAllowedAssignments(d + [c], allowed)
    terms.append(c)
for u, v in B.PHYS: edge(u, v, 0)
for a, b in B.SESS: edge(a, b, I)
md.Add(sum(terms) >= 147 * m)          # flow bound (stronger-rank-screens Sec.2; rank_F2 is a metric)
md.Minimize(sum(terms))
for v in range(B.NV):
    for k in range(nb): md.AddHint(bits[v][k], 0)
sv = cp_model.CpSolver(); sv.parameters.num_workers = 1; sv.parameters.max_time_in_seconds = TL
t0 = time.time(); st = sv.Solve(md)
print(f"m={m} var={var} status={sv.StatusName(st)} best={sv.ObjectiveValue():.0f} bound={sv.BestObjectiveBound():.0f} "
      f"Wm={157*m} flow={147*m} t={time.time()-t0:.0f}s", flush=True)
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    beta = [sum(sv.Value(bits[v][k]) << k for k in range(nb)) for v in range(B.NV)]
    print("check s =", B.s_of(beta, tab, I), " distinct frames used:", len(set(beta)))
