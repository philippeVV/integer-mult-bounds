"""Singer cyclic orders on full PG(2,q): P_i at 2(m i mod N), L_l at 2((m l - c) mod N)+1, circle cut before point 0.
cost(m,c) = sum_x (k - #{y >= N-x})^+ + sum_z (k - #{y <= z})^+,  y_d = (m d - c) mod N.
Scans all m coprime to N and all c; writes best order(s) and checks a few cuts explicitly.  usage: cyc.py q"""
import sys, numpy as np
from math import gcd
q=int(sys.argv[1]); D=np.array(list(map(int,open(f'D{q}.txt').read().split()))); N=q*q+q+1; k=(q+1)//2
Sval=((q+1)*(q*q+q+1)-int(round(q**.5))*(q*q+1))//4
res=[]
cs=np.arange(N)
for m in range(1,N):
    if gcd(m,N)!=1: continue
    Y=(m*D[None,:]-cs[:,None])%N             # (c, d)
    H=np.zeros((N,N),dtype=np.int32)
    np.add.at(H,(np.repeat(cs,len(D)),Y.ravel()),1)
    le=np.cumsum(H,axis=1)                    # #{y<=z}
    tot=len(D)
    ge=tot-np.concatenate([np.zeros((N,1),dtype=np.int32),le[:,:-1]],axis=1)  # #{y>=t}
    # point x: count #{y >= N-x}; x=0 -> y>=N -> 0
    jP=np.concatenate([np.zeros((N,1),dtype=np.int32),ge[:,::-1][:,:-1]],axis=1)  # x=1..N-1 -> ge[N-x]
    cost=np.maximum(k-jP,0).sum(1)+np.maximum(k-le,0).sum(1)
    c=int(cost.argmin()); res.append((int(cost[c]),m,c))
res.sort(); print('q',q,'S',Sval,'best cyclic',res[:8])
def order(m,c):
    keys=[(2*((m*i)%N),i) for i in range(N)]+[(2*((m*l-c)%N)+1,N+l) for l in range(N)]
    return [v for _,v in sorted(keys)]
nb=[list(map(int,x.split())) for x in open(f'graph{q}.txt').read().split('\n')[1:2*N+1]]
def cost_of(o):
    pos=[0]*(2*N)
    for t,v in enumerate(o): pos[v]=t
    return sum(max(0,k-sum(pos[u]<pos[v] for u in nb[v])) for v in range(2*N))
for cst,m,c in res[:3]:
    o=order(m,c); cc=[cost_of(o[t:]+o[:t]) for t in range(0,2*N,max(1,2*N//40))]
    print('m',m,'c',c,'formula',cst,'explicit cut0',cc[0],'min over sampled cuts',min(cc))
best=res[0]; txt='\n'.join(map(str,order(best[1],best[2])))+'\n'
import os
if os.path.exists(f'cyc{q}.txt'): print(f'cyc{q}.txt identical to this run:', open(f'cyc{q}.txt').read()==txt)
else: open(f'cyc{q}.txt','w').write(txt)
F={9:166,25:3480,49:26591}[q]   # expXQ: best cyclic F = cost (class-sign orientation is Eulerian)
print(f'best cyclic cost {best[0]} = S + {best[0]-Sval} (expected {F})')
print(('PASS' if best[0]==F else 'FAIL')+f'_EXPXQ_CYC{q}')
