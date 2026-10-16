"""Lexicographic quotient layouts on the full PG(2,q): classes = i mod n' (n' | N), class sequence x_t = s t mod n',
point x+n' t at key (T, 2(m t mod b)), line j = x + r + n' t at key (T, 2((m t - cc) mod b)+lo), b = N/n'.
Scans s, r, m, cc, lo (optionally random subsample).  n'=c gives Baer layouts.  usage: quot.py q n' [maxcombos]"""
import sys, math, numpy as np, random
q=int(sys.argv[1]); n1=int(sys.argv[2]); MAXC=int(sys.argv[3]) if len(sys.argv)>3 else 10**9
r_=math.isqrt(q); N=q*q+q+1; k=(q+1)//2; b=N//n1; assert n1*b==N
L=open(f'graph{q}.txt').read().split('\n'); V=2*N; nb=np.array([list(map(int,L[1+v].split())) for v in range(V)])
Sval=((q+1)*(q*q+q+1)-r_*(q*q+1))//4
def bcost(o):
    pos=np.empty(V,dtype=np.int64); pos[o]=np.arange(V)
    j=(pos[nb]<pos[:,None]).sum(1); return int(np.maximum(k-j,0).sum())
SONE=len(sys.argv)>4
combos=[(s_,r,m,cc,lo) for s_ in ([1] if SONE else range(1,n1)) if math.gcd(s_,n1)==1 for r in range(n1) for m in range(1,b) if math.gcd(m,b)==1
        for cc in range(b) for lo in (-1,1)]
if len(combos)>MAXC: random.seed(1); combos=random.sample(combos,MAXC)
t=np.arange(b); best=None; import heapq; top=[]
for s_,r,m,cc,lo in combos:
    Tpos=np.empty(n1,dtype=np.int64); Tpos[(s_*np.arange(n1))%n1]=np.arange(n1)
    key=np.empty(V,dtype=np.int64)
    x=np.arange(n1)[:,None]
    key[(x+n1*t[None,:]).ravel()]=(Tpos[:,None]*10*b+2*((m*t[None,:])%b)).ravel()
    key[(N+(x+r+n1*t[None,:])%N).ravel()]=(Tpos[:,None]*10*b+2*((m*t[None,:]-cc)%b)+lo).ravel()
    o=np.argsort(key,kind='stable'); cst=bcost(o)
    heapq.heappush(top,(-cst,random.random(),o.copy()))
    if len(top)>40: heapq.heappop(top)
    if best is None or cst<best[0]: best=(cst,(s_,r,m,cc,lo),o); print('best',cst,best[1],flush=True)
print(f"q {q} n' {n1} S {Sval} best quotient layout {best[0]} params {best[1]} (combos {len(combos)})")
open(f'quot{q}_{n1}.txt','w').write('\n'.join(map(str,best[2]))+'\n')

import importlib.util
sys.argv=[sys.argv[0],f'graph{q}.txt']
spec=importlib.util.spec_from_file_location('mcf','mcf.py'); M=importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
bestF=None
for negc,_,o in top:
    c_,Fv=M.F(o)
    if bestF is None or Fv<bestF[0]: bestF=(Fv,c_,o)
print(f"best F (min back arcs) among top-{len(top)} cost layouts: F={bestF[0]} (cost {bestF[1]})")
open(f'quotF{q}_{n1}.txt','w').write('\n'.join(map(str,bestF[2]))+'\n')
