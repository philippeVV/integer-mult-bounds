"""Hadamard-class family: N = (+)_l (G/2)[C_l] (frame F0 terms per class) + private terms (cost 2) for class 0.
Constraint: no lost pair (odd, not plane-orth) inside a class l>=1. Cost = 2|C0| + sum_l sum_k formrank(planes(C_l & D[k])).
SA over r and classes. Exact final check of the decoder by verify()."""
import sys, json, random, math, os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','lower'))
from base import *
from fractions import Fraction as Fr
emap={int(k):v for k,v in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','lower','e8map.json'))).items()}
ip=lambda a,b: sum(x*y for x,y in zip(a,b))//4
Rt=[emap[q] for q in P]; frame=[emap[x] for x in sorted(set(emap)-set(P))]
c=[[ip(Rt[a],f) for f in frame] for a in range(n)]
D=[[a for a in range(n) if c[a][k]] for k in range(8)]
Ks=[[k for k in range(8) if c[a][k]] for a in range(n)]
ODD=[[b for b in range(n) if b!=a and dot(P[a],P[b])==0] for a in range(n)]
RC=[rchoices(q) for q in P]
def lost(r,a,b): return dot(P[a],r[b]) or dot(r[a],P[b]) or dot(r[a],r[b])
def frank(vs): return form_rank(vs)
def run(seed, ncls, iters, T0=2.0, T1=0.05, lam=3.0):
    rng=random.Random(seed)
    r=[rng.choice(RC[a]) for a in range(n)]
    cl=[rng.randrange(ncls+1) for a in range(n)]
    def termcost(l,k):
        vs=[w for a in D[k] if cl[a]==l for w in (P[a],r[a])]
        return frank(vs) if vs else 0
    TC={(l,k):termcost(l,k) for l in range(1,ncls+1) for k in range(8)}
    def viol_port(a):
        if cl[a]==0: return 0
        return sum(1 for b in ODD[a] if cl[b]==cl[a] and lost(r,a,b))
    def total():
        v=sum(viol_port(a) for a in range(n))//2
        return 2*sum(1 for a in range(n) if cl[a]==0)+sum(TC.values()), v
    cost,v=total(); E=cost+lam*v; best=None
    for it in range(iters):
        T=T0*(T1/T0)**(it/iters)
        a=rng.randrange(n); old=(r[a],cl[a])
        v0=viol_port(a); affected=set()
        if rng.random()<0.5: r[a]=rng.choice(RC[a])
        else: cl[a]=rng.randrange(ncls+1)
        for l in {old[1],cl[a]}:
            if l>0:
                for k in Ks[a]: affected.add((l,k))
        oldTC={x:TC[x] for x in affected}
        for x in affected: TC[x]=termcost(*x)
        v1=viol_port(a)
        dcost=sum(TC[x]-oldTC[x] for x in affected)+2*((cl[a]==0)-(old[1]==0))
        # violations involving a change by v1-v0 (pairs counted once)
        # but changing r[a] doesn't affect other pairs; class change only a's pairs: ok
        dE=dcost+lam*(v1-v0)
        if dE<=0 or rng.random()<math.exp(-dE/T):
            cost+=dcost; v+=v1-v0; E+=dE
            if v==0 and (best is None or cost<best[0]): best=(cost,list(r),list(cl))
        else:
            r[a],cl[a]=old
            for x in affected: TC[x]=oldTC[x]
    return best
def build_terms(r,cl,ncls):
    terms=[]
    for a in range(n):
        if cl[a]==0: terms.append(({a:Fr(1)},{a:Fr(1)}))
    for l in range(1,ncls+1):
        for k in range(8):
            S=[a for a in D[k] if cl[a]==l]
            if S: terms.append(({a:Fr(c[a][k],4) for a in S},{a:Fr(c[a][k]) for a in S}))
    return terms
if __name__=='__main__':
    seed=int(sys.argv[1]); ncls=int(sys.argv[2]); iters=int(sys.argv[3])
    b=run(seed,ncls,iters)
    if b:
        bad,C,costs=verify(b[1],build_terms(b[1],b[2],ncls))
        print('seed',seed,'ncls',ncls,'SA cost',b[0],'exact verify bad',bad,'C',C,'class sizes',[b[2].count(l) for l in range(ncls+1)],flush=True)
        json.dump({'r':b[1],'cl':b[2],'C':C,'bad':bad},open('cls_%d_%d.json'%(seed,ncls),'w'))
    else: print('seed',seed,'no feasible',flush=True)
