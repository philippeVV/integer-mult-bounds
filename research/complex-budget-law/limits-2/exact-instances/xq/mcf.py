"""F(sigma) = min #back arcs over all Eulerian (k-in/k-out) orientations, for a given order (ortools min-cost flow).
F >= Boot cost(sigma); F = cost iff the back-arc transportation is feasible.  F(sigma) is the session count of the best
helper-free schedule with topological order sigma.  usage: mcf.py graph.txt order.txt [order2 ...]"""
import sys, numpy as np
from ortools.graph.python import min_cost_flow
L=open(sys.argv[1]).read().split('\n'); V,d=map(int,L[0].split()); k=d//2
nb=np.array([list(map(int,L[1+v].split())) for v in range(V)])
ed=np.array([(a,b) for a in range(V) for b in nb[a] if a<b])
def F(order):
    pos=np.empty(V,dtype=np.int64); pos[np.asarray(order)]=np.arange(V)
    a,b=ed[:,0],ed[:,1]; fw=pos[a]<pos[b]; u=np.where(fw,a,b); w=np.where(fw,b,a)
    j=np.bincount(w,minlength=V); cost=int(np.maximum(k-j,0).sum())
    # start: all edges forward (u->w). flipping edge e (w->u) costs 1 back arc; in-degree target k
    m=min_cost_flow.SimpleMinCostFlow()
    # flow along edge w->u means flip: decreases indeg w, increases indeg u
    m.add_arcs_with_capacity_and_unit_cost(w,u,np.ones(len(u),np.int64),np.ones(len(u),np.int64))
    sup=j-k   # vertices with excess indegree must push out flips
    for v in range(V): m.set_node_supply(v,int(sup[v]))
    assert m.solve()==m.OPTIMAL
    return cost,int(m.optimal_cost())
for f in (sys.argv[2:] if __name__=="__main__" else []):
    o=[int(x) for x in open(f).read().split()]; c,Fv=F(o); print(f"{f}: Boot cost {c}  min back arcs F = {Fv}  (excess over cost {Fv-c})")
