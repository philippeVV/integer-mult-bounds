"""Exact causal schedule search on H(q) incidence (Eulerian K-in/K-out orientation + order, min |R|).
Start: Fiedler order (and BFS-ball order); SA over order with exact min-cost-flow orientation (expNI mcf_eval)."""
import sys, random, math, pickle, numpy as np
from ortools.graph.python import min_cost_flow
q=int(sys.argv[1]); seed=int(sys.argv[2]); iters=int(sys.argv[3]); T0=float(sys.argv[4])
A=np.load(f"A{q}.npy"); V=A.shape[0]; K=int(A[0].sum())//2
nbrs=[list(np.nonzero(A[v])[0]) for v in range(V)]; edges=[(a,b) for a in range(V) for b in nbrs[a] if a<b]
def mcf_eval(pos,want=False,BIG=10000):
    mcf=min_cost_flow.SimpleMinCostFlow(); hub=V; arcs=[]
    for (x,y) in edges:
        u,w=(x,y) if pos[x]<pos[y] else (y,x)
        arcs.append(mcf.add_arc_with_capacity_and_unit_cost(int(u),int(w),1,1))
    for v in range(V):
        b=K-sum(1 for u in nbrs[v] if pos[u]<pos[v]); mcf.set_node_supply(v,b)
        if b>0: mcf.add_arc_with_capacity_and_unit_cost(v,hub,b,BIG)
        elif b<0: mcf.add_arc_with_capacity_and_unit_cost(hub,v,-b,BIG)
    assert mcf.solve()==mcf.OPTIMAL
    cost=mcf.optimal_cost(); slack=cost//(2*BIG); R=cost-slack*2*BIG
    if want: return slack,R,[ei for ei,a in enumerate(arcs) if mcf.flow(a)>0]
    return slack,R
rng=random.Random(seed)
L=np.diag(A.sum(1).astype(float))-A.astype(float); w,U=np.linalg.eigh(L)
f=U[:,1]+(U[:,2]*rng.uniform(-.3,.3) if seed else 0)
order=list(np.argsort(f))
def ev(order):
    p=[0]*V
    for t,v in enumerate(order): p[v]=t
    s,R=mcf_eval(p); return 10000*s+R
cur=ev(order); best=cur; bo=order[:]; print("start",cur,flush=True)
for it in range(iters):
    T=T0*(1-it/iters)+1e-3
    i=rng.randrange(V); Lj=rng.choice([1,2,3,5,8,13,30,60])
    j=min(V-1,max(0,i+rng.randint(-Lj,Lj)))
    if i==j: continue
    new=order[:]; v=new.pop(i); new.insert(j,v)
    c2=ev(new)
    if c2<=cur or rng.random()<math.exp((cur-c2)/T):
        order=new; cur=c2
        if c2<best: best=c2; bo=order[:]
    if it%5000==0: print(it,cur,best,flush=True)
p=[0]*V
for t,v in enumerate(bo): p[v]=t
s,R,back=mcf_eval(p,want=True)
pickle.dump(dict(pos=p,best=best,back=back,edges=edges),open(f"sched{q}_{seed}.pkl","wb"))
print("final",best,"slack",s,"R",R,flush=True)
