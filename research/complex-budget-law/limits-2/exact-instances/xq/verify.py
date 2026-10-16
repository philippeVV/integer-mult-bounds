"""Gate-by-gate verification of a helper-free causal schedule for the odd-q Hermitian Tanner code on PG(2,q) (expQE verify_sched style).
usage: verify.py p order.txt [dense]
 1. order sigma -> j_v; Boot cost; back-arc transportation (max flow) -> explicit Eulerian orientation (in = out = k at every vertex),
    R = back arcs; asserts Gamma - R is acyclic (all non-R arcs go forward in sigma).
 2. simulate over F_q: session symbols on R arcs, each vertex applies its local bijection (in-arcs -> out-arcs) of its MDS gate
    (point: x_d = G_d a; line: x_d = c(i,d) H_d b, c = theta* on wrap edges -- exactly odd_singer.py's dense model).
    Phi: session symbols -> values recomputed at R-arc tails.  Exact & sound <=> Phi = I <=> dim W = |R|.
 3. optional dense dim W (q=9) for an independent check."""
import sys, os, json, numpy as np
from collections import deque
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)   # python -I omits the script dir; 'code' is ./code.py
sys.argv=[sys.argv[0]]+sys.argv[1:]
import code as C
p,q,N,k=C.p,C.q,C.N,C.k; D=C.D; d=q+1
ordf=os.path.join(HERE, sys.argv[2]); DENSE=len(sys.argv)>3
# ---- F_q as indices 0..q-1 with tables
FQ=sorted(C.Fq); IX={a:i for i,a in enumerate(FQ)}
ADDT=np.array([[IX[int(C.add(a,b))] for b in FQ] for a in FQ],dtype=np.uint8); MULT=np.array([[IX[int(C.mul(a,b))] for b in FQ] for a in FQ],dtype=np.uint8)
NEG=np.array([IX[int(C.neg(a))] for a in FQ],dtype=np.uint8); INV=np.zeros(q,dtype=np.uint8)
for i,a in enumerate(FQ):
    if a: INV[i]=IX[C.inv(a)]
Z,O=IX[0],IX[1]
def mat_inv(A):
    n=A.shape[0]; M=np.concatenate([A,np.full((n,n),Z,dtype=np.uint8)],1)
    for i in range(n): M[i,n+i]=O
    for c in range(n):
        nz=[r for r in range(c,n) if M[r,c]!=Z]
        if not nz: return None
        r=nz[0]; M[[c,r]]=M[[r,c]]; M[c]=MULT[INV[M[c,c]],M[c]]
        for r in range(n):
            if r!=c and M[r,c]!=Z: M[r]=ADDT[M[r],MULT[NEG[M[r,c]],M[c]]]
    return M[:,n:]
def mat_mul(A,B):
    out=np.full((A.shape[0],B.shape[1]),Z,dtype=np.uint8)
    for t in range(A.shape[1]): out=ADDT[out,MULT[A[:,t][:,None],B[t][None,:]]]
    return out
def rank_idx(A):
    M=A.copy(); r=0; R_,Cc=M.shape
    for c in range(Cc):
        if r>=R_: break
        nz=np.nonzero(M[r:,c]!=Z)[0]
        if len(nz)==0: continue
        i=r+nz[0]; M[[r,i]]=M[[i,r]]; M[r]=MULT[INV[M[r,c]],M[r]]
        rows=np.nonzero(M[:,c]!=Z)[0]; rows=rows[rows!=r]
        if len(rows): M[rows]=ADDT[M[rows],MULT[NEG[M[rows,c]][:,None],M[r][None,:]]]
        r+=1
    return r
th=C.gp(C.theta_star_exp); assert th in C.Fq
Gi=np.vectorize(lambda a:IX[int(a)])(C.G); Hi=np.vectorize(lambda a:IX[int(a)])(C.H); thi=IX[th]
# ---- graph: edge e=(i,di): point i (vertex i), line l=(i+d)%N (vertex N+l)
E=N*d; ends=[]; coef=[]
inc=[[] for _ in range(2*N)]; row=[[] for _ in range(2*N)]   # row[v][t]: local generator row (k-vector) of v's t-th edge
for i in range(N):
    for di,dd in enumerate(D):
        l=(i+dd)%N; e=len(ends); ends.append((i,N+l)); inc[i].append(e); inc[N+l].append(e)
        row[i].append(Gi[di]); c=thi if i+dd>=N else O
        row[N+l].append(MULT[c,Hi[di]])
Mloc=[np.array(r,dtype=np.uint8) for r in row]
import os
if os.environ.get('PERTURB'):   # negative control: rescale one coordinate of line codes (breaks the global code, keeps MDS gates)
    mode=os.environ['PERTURB']
    for v in range(N,2*N):
        if mode=='one' and v!=N: continue
        Mloc[v][0]=MULT[IX[C.gp(N)],Mloc[v][0]]          # g^N generates F_q^*

order=[int(x) for x in open(ordf).read().split()]; V=2*N; assert sorted(order)==list(range(V))
pos=[0]*V
for t,v in enumerate(order): pos[v]=t
other=lambda e,v: ends[e][1] if ends[e][0]==v else ends[e][0]
j=[sum(pos[other(e,v)]<pos[v] for e in inc[v]) for v in range(V)]
cost=sum(max(0,k-x) for x in j); print(f"q={q} V={V} E={E} k={k}  Boot cost of order = {cost}",flush=True)
# ---- transportation: surplus v (j>k) sends j-k back arcs to earlier deficit nbrs u (j<k), u receives k-j
G={}; cap={}
def addE(a,b,c):
    G.setdefault(a,[]).append(b); G.setdefault(b,[]).append(a); cap[(a,b)]=cap.get((a,b),0)+c; cap.setdefault((b,a),0)
need=0
for v in range(V):
    if j[v]>k:
        addE('s',('o',v),j[v]-k); need+=j[v]-k
        for e in inc[v]:
            u=other(e,v)
            if pos[u]<pos[v] and j[u]<k: addE(('o',v),('e',e),1); addE(('e',e),('i',u),1)
    if j[v]<k: addE(('i',v),'t',k-j[v])
flow=0
while True:
    par={'s':None}; dq=deque(['s'])
    while dq and 't' not in par:
        a=dq.popleft()
        for b in G.get(a,[]):
            if b not in par and cap[(a,b)]>0: par[b]=a; dq.append(b)
    if 't' not in par: break
    b='t'
    while par[b] is not None: cap[(par[b],b)]-=1; cap[(b,par[b])]+=1; b=par[b]
    flow+=1
print(f"back-arc transportation flow {flow}/{need}",flush=True)
if flow!=need: print("NOT SCHEDULABLE FROM THIS ORDER (transportation infeasible)"); sys.exit(1)
back=set(e for e in range(E) if ('e',e) in G and cap[(('o',None),)] if False) if False else set()
for v in range(V):
    for e in inc[v]:
        if (('o',v),('e',e)) in cap and cap[(('o',v),('e',e))]==0 and pos[other(e,v)]<pos[v]: back.add(e)
arcs=[]
for e,(a,b) in enumerate(ends):
    lo,hi=(a,b) if pos[a]<pos[b] else (b,a)
    arcs.append((hi,lo) if e in back else (lo,hi))
indeg=[0]*V; outdeg=[0]*V
for a,b in arcs: outdeg[a]+=1; indeg[b]+=1
assert all(x==k for x in indeg) and all(x==k for x in outdeg), "not Eulerian k/k"
R=sorted(back); assert len(R)==cost
assert all(pos[a]<pos[b] for e,(a,b) in enumerate(arcs) if e not in back)   # Gamma-R acyclic, sigma is a topological order
print(f"Eulerian {k}-in/{k}-out orientation OK; |R| = {len(R)} back arcs; Gamma-R acyclic (sigma topological)",flush=True)
# ---- simulation
Sn=len(R); Rk={e:t for t,e in enumerate(R)}
val=np.full((E,Sn),Z,dtype=np.uint8); have=[False]*E
for e,t in Rk.items(): val[e,t]=O; have[e]=True
Phi=np.full((Sn,Sn),Z,dtype=np.uint8); bad=0; seen_types=set()
for v in order:
    ins=[t for t,e in enumerate(inc[v]) if arcs[e][1]==v]; outs=[t for t,e in enumerate(inc[v]) if arcs[e][0]==v]
    assert len(ins)==k==len(outs)
    Ai=mat_inv(Mloc[v][ins])
    if Ai is None: bad+=1; continue
    T=mat_mul(Mloc[v][outs],Ai)            # x_O = M_O M_I^-1 x_I
    ein=[inc[v][t] for t in ins]; assert all(have[e] for e in ein)
    out=mat_mul(T,val[ein])
    for r,t in enumerate(outs):
        e=inc[v][t]
        if e in Rk: Phi[:,Rk[e]]=out[r]
        else: val[e]=out[r]; have[e]=True
print(f"gates: {V} simulated, non-invertible: {bad}",flush=True)
Dm=Phi.T.copy()
for t in range(Sn): Dm[t,t]=ADDT[Dm[t,t],NEG[O]]
rk=rank_idx(Dm); print(f"rank(Phi - I) = {rk}  => dim W(schedule) = {Sn-rk};  EXACT & SOUND: {rk==0 and bad==0}",flush=True)
if rk==0 and bad==0:   # limits-2: compare with the stored artifact instead of rewriting it (#2838's copy rewrites schedule{q}.json)
    new=json.dumps(dict(q=q,D=D,order=order,arcs=arcs,R=R)); fn=os.path.join(HERE,f'schedule{q}.json')
    same=os.path.exists(fn) and open(fn).read()==new
    print(f"schedule{q}.json byte-identical to this run: {same}",flush=True)
    if not os.path.exists(fn): open(fn,'w').write(new)
if DENSE:
    Mx=np.full((E,2*N*k),Z,dtype=np.uint8)
    for v in range(V):
        for t,e in enumerate(inc[v]):
            if v<N: Mx[e,v*k:(v+1)*k]=Mloc[v][t]
            else: Mx[e,v*k:(v+1)*k]=NEG[Mloc[v][t]]
    dW=2*N*k-rank_idx(Mx); print("dense dim W =",dW)
EXPECT={'':({9:0,25:31,49:157}.get(q),None),'one':(1,165),'all':(61,105)}[os.environ.get('PERTURB','')]   # cyclic excess; negative controls (q = 9)
okk = rk==EXPECT[0] and bad==0 and (not DENSE or EXPECT[1] is None or dW==EXPECT[1]) and (not DENSE or EXPECT[1] is not None or dW==Sn-EXPECT[0])
print(('PASS' if okk else 'FAIL')+f'_EXPXQ_VERIFY{q}'+('_'+os.environ['PERTURB'].upper() if os.environ.get('PERTURB') else ''))
