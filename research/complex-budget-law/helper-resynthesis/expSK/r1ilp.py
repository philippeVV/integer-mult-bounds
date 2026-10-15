"""Per-level lower bound on helper log-mass under R1 (targets absorb only), edge form, exact MILP (HiGHS).
Model (favourable relaxation of the literal audit; OR-reach, time-free, frames as subspaces of F_2^9):
 * non-target (nt) registers = 112 sources (free, one edge per dim-step 1..7) + helpers (interval of steps a..b, cost f(b-a+1)).
 * E_d = #nt edges at step d (dim d -> d+1) <= 112 + #helpers covering step d   [targets never relay: R1].
 * every pair (s,t) / (s,q) has a delivery node; pairs not absorbed by a target at dim < d pass an nt node at dim d:
   P_d >= 6272 - 112*Mtop(d-1)  (target t absorbs only ports of its own dim-(d-1) frame, <= Mtop(d-1)).
 * node X at dim d with k nt in-edges and j nt out-edges and z absorbing targets carries <= c(k)*(r(j)+z) pairs,
   c(k) = ports covered by k hyperplanes of X (C1, C2, then +C1 each), r(j) = V covered by j hyperplanes of X^perp,
   z <= #ports in X^perp, r+z <= #V in X^perp; #nodes of each profile type <= its number of subspaces (prof*.txt).
 * sum over dim-d nodes: in-edges <= E_{d-1}, out-edges <= E_d, absorbing targets <= 112.
 * cuts: every top needs >= 1 nt arrival (52 > Mtop(7)=32), each U_q needs >= 2 (56 > C1=32) -> E_7 >= 128.
argv: [mode]  mode=R1 (default) | relay (targets relay too: 224 free relays, expSH-style check)."""
import sys, math, json
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
mode=sys.argv[1] if len(sys.argv)>1 else 'R1'
FREE=224 if mode=='relay' else 112
f=lambda w: w*math.log(73/w)
T={d:[list(map(int,l.split()))[1:] for l in open(f'prof{d}.txt')] for d in range(1,9)}
Mtop={d:max(t[0] for t in T[d] if t[4]>=1) for d in T}
AMIN0=int(next((a.split('=')[1] for a in sys.argv if a.startswith('amin=')),'1'))
Mtop={d:(Mtop[d] if d>=AMIN0 else 0) for d in Mtop}
P={1:6272,**{d:6272-112*Mtop[d-1] for d in range(2,9)}}
V=[]; idx={}
def var(name,lo,hi,integ,cost):
    idx[name]=len(V); V.append((lo,hi,integ,cost)); return idx[name]
iv=[(a,b) for a in range(1,8) for b in range(a,8)]
for a,b in iv: var(('x',a,b),0,np.inf,1,f(b-a+1))
for d in range(1,8): var(('E',d),0,np.inf,0,0)
for d in range(1,9): var(('A',d),0,np.inf,0,0)
SPLIT='nosplit' not in sys.argv
AMIN=int(next((a.split('=')[1] for a in sys.argv if a.startswith('amin=')),'1'))  # targets absorb only at dims >= AMIN
rows=[]   # (coeff dict, lo, hi)
for d in range(1,8): rows.append(({idx[('x',a,b)]:-1 for a,b in iv if a<=d<=b}|{idx[('E',d)]:1},-np.inf,FREE))
rows.append(({idx[('E',7)]:1},128 if mode!='relay' else 0,np.inf))
def cfun(t,k,d):
    Pp,C1,C2=t[0],t[1],t[2]
    if d==1: return Pp if k==0 else None
    if k==0: return 0
    return min(Pp,[0,C1,C2][k] if k<=2 else C2+(k-2)*C1)
def rfun(t,j):
    QV,R1,R2=t[3],t[5],t[6]
    return min(QV,[0,R1,R2][j] if j<=2 else R2+(j-2)*R1)
for d in range(1,9):
    inn={}; out={}; tau={}; pairs={}; outp={}; absp={}
    for ti,t in enumerate(T[d]):
        Pp,C1,C2,QV,Qp,R1,R2,cnt=t
        if Pp==0 or QV==0: continue
        cap={}
        Ks=[0] if d==1 else range(1,9)
        for k in Ks:
            c=cfun(t,k,d)
            if c is None or c==0: continue
            if k>=2 and c==cfun(t,k-1,d): break
            Js=[0] if d==8 else range(0,9)
            for j in Js:
                r=rfun(t,j) if d<8 else QV-Qp      # dim 8: centre normal reach is free (U_q register), targets via z
                if j>=1 and r==rfun(t,j-1) and j>1: break
                zmax=min(Qp,QV-r) if mode!='relay' else min(Qp,QV-r)
                y=var(('y',d,ti,k,j),0,np.inf,0,0); cap[y]=1
                inn[y]=k; out[y]=j; pairs[y]=c*r; outp[y]=c*r
                if d<AMIN: zmax=0
                if zmax>0:
                    z=var(('z',d,ti,k,j),0,np.inf,0,0); tau[z]=1; pairs[z]=c; absp[z]=c
                    rows.append(({z:1,y:-zmax},-np.inf,0))
        if cap: rows.append((cap,-np.inf,cnt))
    if d>=2: rows.append(({**{y:k for y,k in inn.items() if k},idx[('E',d-1)]:-1},-np.inf,0))
    if d<=7: rows.append(({**{y:j for y,j in out.items() if j},idx[('E',d)]:-1},-np.inf,0))
    if mode!='relay': rows.append((dict(tau),-np.inf,112))
    else: rows.append((dict(tau),-np.inf,112))
    rows.append((dict(pairs),P[d],np.inf))
    if SPLIT:
        rows.append(({**absp,idx[('A',d)]:-1},0,np.inf))                       # newly absorbed at d <= absorb capacity
        rows.append(({idx[('A',e)]:1 for e in range(1,d+1)},-np.inf,112*Mtop[d] if d<8 else 5824))
        if d<=7: rows.append(({**outp,**{idx[('A',e)]:1 for e in range(1,d+1)}},6272,np.inf))  # rest leaves via out-edges
        else: rows.append(({**{y:v for y,v in pairs.items()},**{idx[('A',e)]:1 for e in range(1,8)}},6272,np.inf))
A=lil_matrix((len(rows),len(V)))
for i,(co,lo,hi) in enumerate(rows):
    for j,v in co.items(): A[i,j]=v
lo=np.array([r[1] for r in rows]); hi=np.array([r[2] for r in rows])
c=np.array([v[3] for v in V]); integ=np.array([v[2] for v in V])
res=milp(c,constraints=LinearConstraint(A.tocsr(),lo,hi),integrality=integ,bounds=Bounds([v[0] for v in V],[v[1] for v in V]),options=dict(time_limit=600,mip_rel_gap=1e-6))
print(mode,'status',res.status,res.message)
xs={f'[{a},{b}]':round(res.x[idx[('x',a,b)]]) for a,b in iv if res.x[idx[('x',a,b)]]>0.5}
E={d:round(res.x[idx[('E',d)]],2) for d in range(1,8)}
print('Mtop',Mtop,'P_d',P)
print('min helper log-mass %.1f'%res.fun,'helpers',sum(xs.values()),xs)
print('nt edges per step',E)
use={}
for (key,i) in idx.items():
    if key[0] in 'yz' and res.x[i]>1e-6:
        d,ti,k,j=key[1:]; t=T[d][ti]
        use.setdefault(d,[]).append((key[0],round(res.x[i],2),'k',k,'j',j,'Pp,C1,C2,QV,Qp,R1,R2,cnt',t))
json.dump(dict(mode=mode,mass=res.fun,helpers=xs,E=E,P=P,Mtop=Mtop,use={str(d):u for d,u in use.items()}),open(f'r1ilp_{mode}{"" if SPLIT else "_nosplit"}_a{AMIN}.json','w'),indent=0)
if '-v' in sys.argv:
    for d in sorted(use):
        for u in sorted(use[d],key=lambda u:-u[1])[:8]: print(d,u)
