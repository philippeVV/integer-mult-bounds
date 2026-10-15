"""Exact: (i) maxports5 = max #ports in a 5-dim subspace (= ports orthogonal to a 4-dim dual Q), branch&bound;
(ii) R6 = max ports(L)*ports(L^perp) over 6-dim L  (= over 3-dim L' by duality), exhaustive over dual triples;
also R6 restricted to ports(L6^perp)>=4 (needed for m>=2) and max ports(L6) given ports(L6^perp)>=4."""
import json, numpy as np, itertools
P=json.load(open('instA.json'))['ports']; dot=lambda a,b:(a&b).bit_count()&1
isport=np.zeros(512,bool); isport[P]=True
M=np.zeros((512,2),dtype=np.uint64)
for q in range(1,512):
    m=0
    for i,p in enumerate(P):
        if dot(p,q)==0: m|=1<<i
    M[q,0]=m&((1<<64)-1); M[q,1]=m>>64
def pc(x):
    if x.shape[0]==0: return np.zeros(0,int)
    return np.unpackbits(x.view(np.uint8).reshape(x.shape[0],-1),axis=1).sum(1)
idx=np.arange(512)
# (ii) over dual triples (a<b<c independent): L6 = <a,b,c>^perp ; ports(L6)=pc(Ma&Mb&Mc); ports(L6^perp)=#ports in span(a,b,c)
R=0; Rm2=0; P6m2=0; best6=0
for a in range(1,512):
    for b in range(a+1,512):
        cs=idx[b+1:]; cs=cs[cs!=(a^b)]
        if len(cs)==0: continue
        x=M[a]&M[b]&M[cs]; n6=pc(x)
        np_=isport[a].astype(int)+isport[b]+isport[a^b]+isport[cs]+isport[cs^a]+isport[cs^b]+isport[cs^a^b]
        prod=n6*np_; R=max(R,prod.max()); best6=max(best6,n6.max())
        sel=np_>=4
        if sel.any(): Rm2=max(Rm2,prod[sel].max()); P6m2=max(P6m2,n6[sel].max())
print('maxports6',best6,'R6',R,' among L6 with >=4 ports in L6^perp: R',Rm2,' max ports(L6)',P6m2,flush=True)
# (i) maxports5 by branch and bound over dual 4-tuples
best=0
for a in range(1,512):
    for b in range(a+1,512):
        mab=M[a]&M[b]
        if pc(mab[None])[0]<=best: continue
        cs=idx[b+1:]
        if len(cs)==0: continue
        x=mab&M[cs]; n=pc(x)
        for c,cnt,xc in zip(cs[n>best],n[n>best],x[n>best]):
            ds=idx[c+1:]; ds=ds[(ds!=(a^b))&(ds!=(a^c))&(ds!=(b^c))&(ds!=(a^b^c))]
            if len(ds)==0: continue
            y=pc(xc&M[ds]); mx=y.max()
            if mx>best: best=int(mx)
print('maxports5',best)
