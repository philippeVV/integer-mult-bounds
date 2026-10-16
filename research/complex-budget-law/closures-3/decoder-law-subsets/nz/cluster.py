import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
# max induced cluster (P3-free) subgraph of the E8 orthogonality graph on a line set L
import sys
from e8 import *
from ortools.sat.python import cp_model
def maxcluster(Ls, tl=600):
    m = cp_model.CpModel(); x = {a: m.NewBoolVar(f"x{a}") for a in Ls}
    S=set(Ls); n=0
    for b in Ls:
        nb=[a for a in Ls if a!=b and IP[a][b]==0]
        for i,a in enumerate(nb):
            for c in nb[i+1:]:
                if IP[a][c]!=0: m.Add(x[a]+x[b]+x[c]<=2); n+=1
    m.Maximize(sum(x.values()))
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=tl; s.parameters.num_workers=1
    r=s.Solve(m)
    return s.StatusName(r), s.ObjectiveValue(), s.BestObjectiveBound(), [a for a in Ls if s.Value(x[a])]
if __name__=="__main__":
    st,v,bd,sol=maxcluster(list(range(120)), float(sys.argv[1]) if len(sys.argv)>1 else 300)
    print(st,v,bd); print(sol)
    print('PASS_EXPNZ_CLUSTER' if st == 'OPTIMAL' and round(v) == 14 else 'FAIL_EXPNZ_CLUSTER')   # cl(P120) = 14: the L2 floor is useless for large Q
