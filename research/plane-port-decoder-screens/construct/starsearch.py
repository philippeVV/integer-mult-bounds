"""local search over r_t in W_t (frame cost 8 kept): minimise sum of row-star costs; also col-star alternative per port (min of row/col orientation is NOT exact-cover safe, so rows only)."""
from csp import *
import json, math
seed=int(sys.argv[1]); iters=int(sys.argv[2]); random.seed(seed)
Cf=[[ip(L[ports[i]],L[f]) for f in fr] for i in range(n)]
dom=[[r for r in range(1,512) if okplane(U[i],r) and all(d(r,labs[fr[k]])==0 for k in range(8) if Cf[i][k])] for i in range(n)]
R=json.load(open('mf_1_1.json'))['R'] if seed==0 else [random.choice(dom[i]) for i in range(n)]
def starcost(s):
    B=[j for j in odd[s] if not orth(s,R[s],j,R[j])]
    if not B: return 0
    vs=[]
    for j in B: vs+= [U[j],R[j]]
    # cost = rank(E_s x span) = 2 - dim(E_s cap span^perp)
    k=0
    for v in (U[s],R[s],U[s]^R[s]):
        if all(d(v,w)==0 for w in vs): k+=1
    return 2 if k==0 else 1  # k counts nonzero vectors of E_s in span^perp: 1 -> dim1 ; 3 -> dim2 (impossible since broken)
cost=[starcost(s) for s in range(n)]; tot=sum(cost); best=tot; T=0.5
for it in range(iters):
    i=random.randrange(n); old=R[i]; R[i]=random.choice(dom[i])
    aff=[i]+odd[i]; newc={s:starcost(s) for s in aff}
    dl=sum(newc[s]-cost[s] for s in aff)
    if dl<=0 or random.random()<math.exp(-dl/T):
        for s in aff: cost[s]=newc[s]
        tot+=dl
        if tot<best: best=tot; json.dump({'R':R},open(f'ss_{seed}.json','w'))
    else: R[i]=old
    T=max(0.02,0.5*(1-it/iters))
print(seed,'best star total',best,'-> C =',64+best,flush=True)
