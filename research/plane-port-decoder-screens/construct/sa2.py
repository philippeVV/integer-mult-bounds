"""tensor(2-class)+drop search. objective 2*drops + W*violations (W ramps up). exact deltas; best-r moves."""
from csp import *
import math, json
seed=int(sys.argv[1]); restrict=int(sys.argv[2]); iters=int(sys.argv[3]); random.seed(seed)
Cf=[[ip(L[ports[i]],L[f]) for f in fr] for i in range(n)]
dom=[]
for i in range(n):
    ds=[r for r in range(1,512) if okplane(U[i],r)]
    if restrict: ds=[r for r in ds if all(d(r,labs[fr[k]])==0 for k in range(8) if Cf[i][k])]
    dom.append(ds)
cls=[random.randrange(2) for _ in range(n)]; R=[random.choice(dom[i]) for i in range(n)]; dr=[0]*n
def viol(i,c,r): return sum(1 for j in odd[i] if not dr[j] and cls[j]==c and not orth(i,r,j,R[j]))
def bestr(i,c):
    b=None
    for r in dom[i]:
        v=viol(i,c,r)
        if b is None or v<b[0]: b=(v,r)
        if v==0: break
    return b
best=None; T=1.0
for it in range(iters):
    W=0.5+4.0*it/iters
    i=random.randrange(n)
    cur=2*dr[i]+(0 if dr[i] else W*viol(i,cls[i],R[i]))
    u=random.random()
    if u<0.1: c,r,x=cls[i],R[i],1-dr[i]
    elif u<0.4: c=1-cls[i]; v,r=bestr(i,c); x=0
    else: c=cls[i]; r=random.choice(dom[i]); x=0
    new=2*x+(0 if x else W*viol(i,c,r))
    dl=new-cur
    if dl<=0 or random.random()<math.exp(-dl/T):
        cls[i],R[i],dr[i]=c,r,x
    T=max(0.05,1.0*(1-it/iters))
    if it%20000==0 or it==iters-1:
        tv=sum(viol(j,cls[j],R[j]) for j in range(n) if not dr[j])//2
        if tv==0 and (best is None or sum(dr)<best):
            best=sum(dr); json.dump({'cls':cls,'R':R,'drop':dr},open(f'sa2_{seed}_{restrict}.json','w'))
print(seed,restrict,'best feasible drops',best,flush=True)
