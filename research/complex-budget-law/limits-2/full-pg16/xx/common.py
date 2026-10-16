import numpy as np
from singer import singer
N,D,_=singer()
V=2*N
K=np.array([8]*N+[9]*N)
nb=[[] for _ in range(V)]
for i in range(N):
    for d in D: l=N+(i+d)%N; nb[i].append(l); nb[l].append(i)
NB=np.array(nb)
ED=np.array([(i,N+(i+d)%N) for i in range(N) for d in D])
def cyc_order(m,c):
    keys=[(2*((m*i)%N),i) for i in range(N)]+[(2*((m*l-c)%N)+1,N+l) for l in range(N)]
    return [v for _,v in sorted(keys)]
def boot_cost(order):
    pos=np.empty(V,dtype=np.int64); pos[np.asarray(order)]=np.arange(V)
    j=(pos[NB]<pos[:,None]).sum(1)
    return int(np.maximum(K-j,0).sum())
def F(order, ret=False):
    from ortools.graph.python import min_cost_flow
    pos=np.empty(V,dtype=np.int64); pos[np.asarray(order)]=np.arange(V)
    a,b=ED[:,0],ED[:,1]; fw=pos[a]<pos[b]; u=np.where(fw,a,b); w=np.where(fw,b,a)
    j=np.bincount(w,minlength=V)
    m=min_cost_flow.SimpleMinCostFlow()
    arcs=m.add_arcs_with_capacity_and_unit_cost(w,u,np.ones(len(u),np.int64),np.ones(len(u),np.int64))
    sup=j-K
    m.set_nodes_supplies(np.arange(V),sup.astype(np.int64))
    assert m.solve()==m.OPTIMAL
    Fv=int(m.optimal_cost())
    if not ret: return Fv
    fl=np.array([m.flow(x) for x in arcs])
    tail=np.where(fl==1,w,u); head=np.where(fl==1,u,w)   # oriented arcs tail->head
    back=fl==1
    return Fv,tail,head,back
