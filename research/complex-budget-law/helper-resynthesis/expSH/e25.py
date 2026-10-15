"""E25 = max ports(A)*ports(B), A perp B, dim A<=2, dim B<=5 (k=1 / single-register bound at levels 3 and 6).
Also mp5 (global max ports in a 5-dim subspace) as a by-product."""
import json, numpy as np, itertools
P=json.load(open('instA.json'))['ports']; n=112; dot=lambda a,b:(a&b).bit_count()&1
M=np.zeros((512,2),np.uint64)
for q in range(512):
    m=sum(1<<i for i,p in enumerate(P) if dot(p,q)==0); M[q,0]=m&((1<<64)-1); M[q,1]=m>>64
pc=lambda x: np.bitwise_count(x).sum(-1)
best=0; arg=None
for a,b in itertools.combinations(range(n),2):
    m=M[P[a]]&M[P[b]]                      # A^perp ports (7-dim)
    S=np.array([0,P[a],P[b],P[a]^P[b]])
    ys=np.arange(512); okY=~np.isin(ys,S)
    Z=(m[None,None,:]&M[:,None,:]&M[None,:,:])
    cnt=pc(Z)
    zz=np.arange(512)
    bad=np.isin(zz[None,:],S)[0][None,:] | np.isin(zz[None,:]^ys[:,None],S) | ~okY[:,None]
    cnt[bad]=0
    v=int(cnt.max())
    if 2*v>best: best=2*v; arg=(a,b,v)
print('E25 (2-port A) <=',best,arg)
