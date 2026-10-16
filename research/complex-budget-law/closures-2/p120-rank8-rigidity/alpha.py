"""Per-component: alpha(T) (exact via CP-SAT) = certified min-rank lower bound (identity principal minor);
span rank of the component's roots = rank of block G/2 = rational upper bound (G/2 restricted is a valid completion: support = F within block = T)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, sys
from ortools.sat.python import cp_model
from e8 import R, F, G, n
from subsys import SUBS, Tgraph
def comps(T):
    comp=-np.ones(n,int); c=0
    for s in range(n):
        if comp[s]>=0: continue
        st=[s]; comp[s]=c
        while st:
            a=st.pop()
            for b in np.nonzero(T[a])[0]:
                if comp[b]<0: comp[b]=c; st.append(b)
        c+=1
    return [np.nonzero(comp==k)[0] for k in range(c)]
def alpha(T, V):
    m=cp_model.CpModel(); x={v:m.NewBoolVar('x%d'%v) for v in V}
    for i,a in enumerate(V):
        for b in V[i+1:]:
            if T[a,b]: m.AddBoolOr([x[a].Not(),x[b].Not()])
    m.Maximize(sum(x.values())); s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=120; s.parameters.num_workers=2
    st=s.Solve(m); return int(s.ObjectiveValue()), int(s.BestObjectiveBound()), st==cp_model.OPTIMAL
if __name__=='__main__':
    lows = []
    for name,Rs in SUBS.items():
        T=Tgraph(Rs)
        # check: within a component, T == F (so block G/2 is supported on T)
        cs=comps(T); tot_lo=0; tot_hi=0; parts=[]
        for V in cs:
            blockF=F[np.ix_(V,V)]; blockT=T[np.ix_(V,V)]
            sr=np.linalg.matrix_rank(R[V].astype(float)) if (blockF==blockT).all() else -1  # -1: G/2 block not supported on T
            al=alpha(T,list(V))
            parts.append((len(V),al[0],sr)); tot_lo+=al[0]; tot_hi+=sr if sr>0 else 999
        print(name,'d=%.2f'%T.sum(1).mean(),'components (size, alpha, spanrank):',parts,' => min-rank in [%d, %d]'%(tot_lo,tot_hi),flush=True)
        if T.sum(1).mean() <= 30: lows.append(tot_lo)
    print('min certified min-rank lower bound over singles with d <= 30:', min(lows))
    print('PASS_ALPHA' if min(lows) >= 12 else 'FAIL_ALPHA')
