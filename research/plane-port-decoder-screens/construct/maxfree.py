"""local search: choose r_t (optionally in W_t, keeping frame-centre cost 8) maximizing plane-orthogonal odd pairs."""
from csp import *
import json, math
seed=int(sys.argv[1]); restrict=int(sys.argv[2]); random.seed(seed)
Cf=[[ip(L[ports[i]],L[f]) for f in fr] for i in range(n)]
dom=[]
for i in range(n):
    ds=[r for r in range(1,512) if okplane(U[i],r)]
    if restrict: ds=[r for r in ds if all(d(r,labs[fr[k]])==0 for k in range(8) if Cf[i][k])]
    dom.append(ds)
R=[random.choice(dom[i]) for i in range(n)]
def fi(i,r): return sum(1 for j in odd[i] if orth(i,r,j,R[j]))
tot=sum(fi(i,R[i]) for i in range(n))//2; best=tot; T=2.0
for it in range(300000):
    i=random.randrange(n); r=random.choice(dom[i])
    dl=fi(i,r)-fi(i,R[i])
    if dl>=0 or random.random()<math.exp(dl/T):
        R[i]=r; tot+=dl
        if tot>best: best=tot; json.dump({'R':R,'free':tot},open(f'mf_{seed}_{restrict}.json','w'))
    T=max(0.05,T*0.99998)
print(seed,restrict,'best free undirected',best,'of 2912; max possible <=',112*28//2)
