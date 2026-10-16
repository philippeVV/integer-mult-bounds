"""Baer-subplane layouts on the full PG(2,q), q=p^2 (Singer model).  Baer classes: points i mod c, c = q - sqrt q + 1;
Baer lines of class x: lines j = x + r0 mod c (r0 = residue holding sqrt q + 1 elements of D).  Layout: class sequence
x_t = x0 + s t mod c; inside a class, cyclic layout of the Baer subplane PG(2,sqrt q): point x+c t at 2(m t mod beta),
line at 2((m t - cc) mod beta)+lo.  Scans s, m, cc, lo; writes best order.  usage: baer.py q"""
import sys, math, numpy as np
from collections import Counter
q=int(sys.argv[1]); r=math.isqrt(q); N=q*q+q+1; k=(q+1)//2; c=q-r+1; beta=q+r+1; assert c*beta==N
D=list(map(int,open(f'D{q}.txt').read().split()))
cnt=Counter(d%c for d in D); r0=max(cnt,key=cnt.get); assert cnt[r0]==r+1
L=open(f'graph{q}.txt').read().split('\n'); V=2*N; nb=np.array([list(map(int,L[1+v].split())) for v in range(V)])
def bcost(order):
    pos=np.empty(V,dtype=np.int64); pos[np.asarray(order)]=np.arange(V)
    j=(pos[nb]<pos[:,None]).sum(1); return int(np.maximum(k-j,0).sum())
Sval=((q+1)*(q*q+q+1)-r*(q*q+1))//4
res=[]
for s_ in range(1,c):
    if math.gcd(s_,c)!=1: continue
    seq=[(s_*t)%c for t in range(c)]
    for m in range(1,beta):
        if math.gcd(m,beta)!=1: continue
        for cc in range(beta):
            for lo in (-1,1):
                key=np.zeros(V)
                for T,x in enumerate(seq):
                    t=np.arange(beta)
                    key[x+c*t]=T*10*beta+2*((m*t)%beta)
                    key[N+(x+r0+c*t)%N]=T*10*beta+2*((m*t-cc)%beta)+lo
                o=np.argsort(key,kind='stable'); res.append((bcost(o),s_,m,cc,lo,o))
res.sort(key=lambda z:z[0]); print('q',q,'S',Sval,'best Baer layouts',[z[:5] for z in res[:8]],flush=True)
open(f'baer{q}.txt','w').write('\n'.join(map(str,res[0][5]))+'\n')
