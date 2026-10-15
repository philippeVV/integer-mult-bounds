"""Min helper log-mass sum f(b-a), f(w)=w ln(73/w), over helper intervals [a,b] (alive, i.e. edge-providing, at dims a..b)
subject to per-dim counts need[d] (helpers alive at dim d). Exact ILP (CP-SAT, scaled costs)."""
import math, sys, json
from ortools.sat.python import cp_model
need=json.loads(sys.argv[1]) if len(sys.argv)>1 else {"2":112,"3":19,"6":19,"8":8}
need={int(k):v for k,v in need.items()}
f=lambda w: w*math.log(73/w)
iv=[(a,b) for a in range(1,9) for b in range(a+1,9)]
m=cp_model.CpModel(); x={i:m.NewIntVar(0,400,str(i)) for i in iv}
for d,k in need.items(): m.Add(sum(x[(a,b)] for a,b in iv if a<=d<=b)>=k)
m.Minimize(sum(int(1e6*f(b-a))*x[(a,b)] for a,b in iv))
s=cp_model.CpSolver(); s.Solve(m)
sol={f'[{a},{b}]':s.Value(x[(a,b)]) for a,b in iv if s.Value(x[(a,b)])}
print('need',need,'-> min helper logmass %.1f'%(s.ObjectiveValue()/1e6),sol)
json.dump(dict(need=need,sol=sol,mass=s.ObjectiveValue()/1e6),open('helperilp.json','w'))
