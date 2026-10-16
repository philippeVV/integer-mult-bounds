"""SA over orders with the TRUE schedule size F(sigma) (min back arcs over Eulerian orientations, min-cost flow) as objective.
Proposals: insertion moves (windowed / next to a neighbour); cheap Boot-cost delta screens out moves with dcost > DC.
usage: saF.py graph.txt init_order secs T0 T1 W DC seed out"""
import sys, math, random, time, numpy as np
import importlib.util
gfile,init,secs,T0,T1,W,DC,seed,out=sys.argv[1],sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]),int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8]),sys.argv[9]
sys.argv=[sys.argv[0],gfile]
spec=importlib.util.spec_from_file_location('mcf','mcf.py'); M=importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
nb=M.nb; V=M.V; k=M.k; rng=random.Random(seed)
order=[int(x) for x in open(init).read().split()]
pos=[0]*V
for t,v in enumerate(order): pos[v]=t
def jv(v): return sum(pos[u]<pos[v] for u in nb[v])
cost,F=M.F(order); best=F; print('start cost',cost,'F',F,flush=True)
t0=time.time(); it=0; ev=0
while time.time()-t0<secs:
    it+=1; T=T0*(T1/T0)**((time.time()-t0)/secs)
    v=rng.randrange(V); a=pos[v]
    if rng.random()<0.5: b=pos[int(nb[v][rng.randrange(len(nb[v]))])]+rng.choice((-1,0,1))
    else: b=a+rng.randint(-W,W)
    b=max(0,min(V-1,b))
    if b==a: continue
    # cost delta
    d=0; dj=0
    for u in nb[v]:
        p=pos[u]; ju=None
        if a<b and a<p<=b: ju=jv(u); d+=max(0,k-ju+1)-max(0,k-ju); dj+=1
        elif b<a and b<=p<a: ju=jv(u); d+=max(0,k-ju-1)-max(0,k-ju); dj-=1
    j0=jv(v); d+=max(0,k-j0-dj)-max(0,k-j0)
    if d>DC: continue
    new=order[:]; new.pop(a); new.insert(b,v)
    c2,F2=M.F(new); ev+=1
    if F2<=F or rng.random()<math.exp(-(F2-F)/T):
        order=new
        for t in range(min(a,b),max(a,b)+1): pos[order[t]]=t
        F=F2; cost=c2
        if F<best:
            best=F; open(out,'w').write('\n'.join(map(str,order))+'\n'); print(f'{time.time()-t0:.0f}s it {it} ev {ev} best F {F} cost {cost} T {T:.3f}',flush=True)
print('final best F',best,'iters',it,'evals',ev)
