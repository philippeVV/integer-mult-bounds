"""Positive control: all-cut-feasible residual class graphs (k=2) should be F2-realisable (SAT), including the
antiparallel (split-class) case. usage: ctl_realise.py LAYOUT ARCS LIMIT"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys
from ortools.sat.python import cp_model
from rg_core import ClassGraph, realise
from max_nres import model
layout=[int(x) for x in sys.argv[1].split(',')]
arcs=[tuple(int(y) for y in a.split('-')) for a in sys.argv[2].split(';')]
limit=int(sys.argv[3])
levels=[]
for t,m in enumerate(layout): levels+=[t]*m
md,d,c,cells=model(levels,arcs,2,target=1)
sols=[]
class CB(cp_model.CpSolverSolutionCallback):
    def on_solution_callback(self):
        sols.append(({k:self.Value(d[k]) for k in cells if self.Value(d[k])},{k:self.Value(c[k]) for k in cells if self.Value(c[k])}))
        if len(sols)>=limit: self.StopSearch()
s=cp_model.CpSolver(); s.parameters.enumerate_all_solutions=True; s.parameters.num_workers=1; s.Solve(md,CB())
st={}
for dd,cc in sols:
    G=ClassGraph(levels,arcs,dd,cc); r,_=realise(G,2,30); key=(G.nres(),r); st[key]=st.get(key,0)+1
print(layout,arcs,len(sols),'(n_res,F2) ->',st)
print('PASS_EXPRG_CTL' if set(v for _, v in st) == {'SAT'} else 'FAIL_EXPRG_CTL')   # positive control: every all-cut-feasible class graph is F2-realisable
