"""All residual configs with 1 <= n_res <= k passing balance + all-L cuts: are they residually routable?
usage: route_lek.py K LAYOUT"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys
from ortools.sat.python import cp_model
from max_nres import model, placements
from route import residual_routable
k=int(sys.argv[1]); layout=[int(x) for x in sys.argv[2].split(',')]
levels=[]
for t,m in enumerate(layout): levels+=[t]*m
tot=0; bad=[]
for arcs in placements(levels,k):
    md,d,c,cells=model(levels,arcs,2,target=1)
    md.Add(sum(d.values())<=k)
    sols=[]
    class CB(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            sols.append(({kk:self.Value(d[kk]) for kk in cells if self.Value(d[kk])},{kk:self.Value(c[kk]) for kk in cells if self.Value(c[kk])}))
    s=cp_model.CpSolver(); s.parameters.enumerate_all_solutions=True; s.parameters.num_workers=1; s.Solve(md,CB())
    for dd,cc in sols:
        tot+=1
        if not residual_routable(levels,arcs,dd,cc): bad.append((arcs,dd,cc))
print(f"k={k} layout={layout}: {tot} cut-feasible configs with 1<=n_res<=k; residual NOT routable: {len(bad)}",flush=True)
for b in bad[:5]: print('  ',b)
print('PASS_EXPRG_ROUTE' if not bad else 'FAIL_EXPRG_ROUTE')   # every cut-feasible config with n_res <= k routes
