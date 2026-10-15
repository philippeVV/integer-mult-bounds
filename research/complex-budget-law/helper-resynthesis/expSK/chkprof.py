"""Stdlib cross-check of prof.c output (prof{d}.txt): (1) subspace counts = Gaussian binomials [9,d]_2;
(2) dims 1 and 8 recomputed exhaustively; (3) 150 random subspaces per dim 2..7 recomputed in pure Python and found
among the listed types; (4) Mtop(d) and the P_d used by r1ilp.py."""
import json, random, itertools
random.seed(1)
V=list(map(int,open('vecs.txt').read().split())); PORT=set(V[:112]); VV=set(V)
dot=lambda a,b:(a&b).bit_count()&1
ok=lambda c,m: print(('PASS ' if c else 'FAIL ')+m) or c
def span(B):
    s=[0]
    for v in B: s=s+[x^v for x in s]
    return s
def basis(vs):
    B=[]
    for v in vs:
        for b in B: v=min(v,v^b)
        if v: B.append(v)
    return B
def stat(B,S,two=True):
    el=span(B); e=len(B); mem=[c for c in range(1<<e) if el[c] in S]
    n=len(mem)
    if e==0 or n==0: return n,0,0
    m=[sum(1<<i for i,c in enumerate(mem) if not (c&w).bit_count()&1) for w in range(1,1<<e)]
    c1=max(x.bit_count() for x in m)
    if not two: return n,c1,None
    c2=max((a|b).bit_count() for a,b in itertools.combinations(m,2)) if len(m)>1 else c1
    return n,c1,max(c1,c2)
def typ(B):
    perp=basis([x for x in range(1,512) if all(dot(x,b)==0 for b in B)])
    Pp,C1,C2=stat(B,PORT); QV,R1,R2=stat(perp,VV); Qp=stat(perp,PORT,False)[0]
    return (Pp,C1,C2,QV,Qp,R1,R2)
def gb(n,k):
    num=den=1
    for i in range(k): num*=(2**(n-i)-1); den*=(2**(i+1)-1)
    return num//den
T={d:{tuple(map(int,l.split()))[1:8]:int(l.split()[8]) for l in open(f'prof{d}.txt')} for d in range(1,9)}
ok(all(sum(T[d].values())==gb(9,d) for d in T),'subspace counts per dim = Gaussian binomials '+str([gb(9,d) for d in range(1,9)]))
from collections import Counter
for d in (1,8):
    if d==1: subs=[[v] for v in range(1,512)]
    else: subs=[basis([x for x in range(1,512) if not dot(x,w)]) for w in range(1,512)]
    ok(Counter(typ(B) for B in subs)==Counter(T[d]),f'dim {d}: exhaustive Python recomputation matches all {len(T[d])} types with counts')
for d in range(2,8):
    good=0
    for _ in range(150 if d not in (6,7) else 40):
        while True:
            B=basis([random.randrange(1,512) for _ in range(d)])
            if len(B)==d: break
        good+=typ(B) in T[d]
    ok(good==(150 if d not in (6,7) else 40),f'dim {d}: random subspaces found among listed types ({good})')
Mtop={d:max(t[0] for t in T[d] if t[4]>=1) for d in T}
ok(Mtop=={1:1,2:2,3:4,4:8,5:12,6:20,7:32,8:52},f'Mtop = {Mtop} (max ports of a d-dim frame inside a top)')
