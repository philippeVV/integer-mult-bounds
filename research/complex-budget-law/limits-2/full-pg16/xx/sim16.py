"""Gate-by-gate F16 simulation of a causal schedule for the full PG(2,16) Hermitian code (non-square gates: points 8->9, lines 9->8).
usage: sim16.py schedfile [json_out]   schedfile = 'order' line, '#arcs' line, then 'tail head' arcs (seed format).
Checks: arcs = Singer incidences once each; in-degree 8 (points) / 9 (lines); R = arcs backward in the order;
Gamma-R acyclic by Kahn (order-independent); every gate's in-set is an information set of its local code;
Phi: session symbols on R -> values recomputed at R tails; A = Phi - I; dim W(schedule) = |R| - rank A.
Completion: dead sessions C must cover the support of A (A[S,S]=0 for S=R\\C): min cover via CP-SAT."""
import sys, os, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)   # python -I omits the script dir
from code16 import make
C=make(); N=C['N']; nb=C['nb']; rows=C['rows']; ADD=C['add']; MUL=C['mul']; INV=C['inv']; V=2*N
K=[8]*N+[9]*N
L=open(os.path.join(HERE, sys.argv[1])).read().split('\n'); order=[int(x) for x in L[0].split()]; ne=int(L[1])
arcs=[tuple(map(int,L[2+t].split())) for t in range(ne)]
assert sorted(order)==list(range(V))
inc=sorted((i,N+(i+d)%N) for i in range(N) for d in C['D'])
assert sorted(tuple(sorted(a)) for a in arcs)==inc, "arcs are not the incidences"
indeg=[0]*V
for a,b in arcs: indeg[b]+=1
assert indeg==K, "in-degrees not 8/9"
pos=[0]*V
for t,v in enumerate(order): pos[v]=t
R=[e for e,(a,b) in enumerate(arcs) if pos[a]>pos[b]]; Rk={e:t for t,e in enumerate(R)}; Sn=len(R)
# Kahn on Gamma-R
adj=[[] for _ in range(V)]; dg=[0]*V
for e,(a,b) in enumerate(arcs):
    if e not in Rk: adj[a].append(b); dg[b]+=1
st=[v for v in range(V) if dg[v]==0]; seen=0
while st:
    v=st.pop(); seen+=1
    for w in adj[v]:
        dg[w]-=1
        if dg[w]==0: st.append(w)
print(f"structure: incidences OK, in-degrees 8/9 OK, |R|={Sn}, Gamma-R acyclic (Kahn): {seen==V}",flush=True)
assert seen==V
def minv(A):
    n=A.shape[0]; M=np.concatenate([A,np.eye(n,dtype=np.uint8)],1)
    for c in range(n):
        nz=[r for r in range(c,n) if M[r,c]]
        if not nz: return None
        r=nz[0]; M[[c,r]]=M[[r,c]]; M[c]=MUL[INV[M[c,c]],M[c]]
        for r in range(n):
            if r!=c and M[r,c]: M[r]=ADD[M[r],MUL[M[r,c],M[c]]]
    return M[:,n:]
def mmul(A,B):
    out=np.zeros((A.shape[0],B.shape[1]),dtype=np.uint8)
    for t in range(A.shape[1]): out=ADD[out,MUL[A[:,t][:,None],B[t][None,:]]]
    return out
def rank(A):
    M=A.copy(); r=0; nr,nc=M.shape
    for c in range(nc):
        if r>=nr: break
        nz=np.nonzero(M[r:,c])[0]
        if len(nz)==0: continue
        i=r+nz[0]; M[[r,i]]=M[[i,r]]; M[r]=MUL[INV[M[r,c]],M[r]]
        rr=np.nonzero(M[:,c])[0]; rr=rr[rr!=r]
        if len(rr): M[rr]=ADD[M[rr],MUL[M[rr,c][:,None],M[r][None,:]]]
        r+=1
    return r
eid={}
for e,(a,b) in enumerate(arcs): eid[(a,b)]=e; eid[(b,a)]=e
val={e:None for e in range(ne)}
for e,t in Rk.items(): x=np.zeros(Sn,dtype=np.uint8); x[t]=1; val[e]=x
Phi=np.zeros((Sn,Sn),dtype=np.uint8); bad=0
for v in order:
    ins=[t for t,u in enumerate(nb[v]) if arcs[eid[(v,u)]][1]==v]; outs=[t for t,u in enumerate(nb[v]) if arcs[eid[(v,u)]][0]==v]
    assert len(ins)==K[v] and len(outs)==17-K[v]
    Ai=minv(rows[v][ins])
    if Ai is None: bad+=1; continue
    T=mmul(rows[v][outs],Ai)
    X=np.array([val[eid[(v,nb[v][t])]] for t in ins])
    assert all(x is not None for x in X)
    Y=mmul(T,X)
    for r,t in enumerate(outs):
        e=eid[(v,nb[v][t])]
        if e in Rk: Phi[:,Rk[e]]=Y[r]
        else: val[e]=Y[r]
print(f"gates simulated: {V} (points 8->9, lines 9->8), singular gate in-sets: {bad}",flush=True)
A=Phi.T.copy()                       # A[t_out, s_in]
for t in range(Sn): A[t,t]^=1        # char 2: Phi - I = Phi + I
rk=rank(A); print(f"rank(Phi - I) = {rk}  => dim W(schedule) = {Sn-rk}  (code dim 900 expected)",flush=True)
nzr=np.nonzero(A.any(1))[0]; nzc=np.nonzero(A.any(0))[0]
print(f"support of Phi-I: {len(nzr)} nonzero rows (outputs), {len(nzc)} nonzero cols (inputs), nnz={int((A!=0).sum())}",flush=True)
cover=None
if rk>0:
    from ortools.sat.python import cp_model
    m=cp_model.CpModel(); x=[m.NewBoolVar(f'x{t}') for t in range(Sn)]
    ii,jj=np.nonzero(A)
    for i,j in set(zip(ii.tolist(),jj.tolist())):
        if i==j: m.Add(x[i]==1)
        else: m.AddBoolOr([x[i],x[j]])
    m.Minimize(sum(x)); s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=120; s.parameters.num_workers=1
    st_=s.Solve(m); cover=[R[t] for t in range(Sn) if s.Value(x[t])]
    c=len(cover); S=Sn-c
    print(f"min dead-session cover c = {c} ({s.StatusName(st_)}): exact sessions {S}, edges {ne}-{c}, s = {S}/{ne-c} = {S/(ne-c):.5f}; "
          f"R4 formula (|R|-dim) -> c={rk}: s_c=(900-{rk})/({ne}-2*{rk}) = {(900-rk)/(ne-2*rk):.5f}",flush=True)
    # re-simulate the cut circuit: sessions R\C, C zero-born at heads, erased at tails -> Phi[S,S] must be I
    keep=[t for t in range(Sn) if R[t] not in set(cover)]
    sub=A[np.ix_(keep,keep)]; print(f"check (Phi - I)[S,S] == 0: {not sub.any()}",flush=True)
if len(sys.argv)>2:   # limits-2: compare with the stored artifact (cover may differ between optimal solutions; compare the rest)
    J=json.load(open(os.path.join(HERE, sys.argv[2])))
    same=all(J[k]==v for k,v in dict(D=C['D'],order=order,arcs=[list(a) for a in arcs],R=R,rank_PhiI=rk).items()) and len(J['dead'])==len(cover)
    print(f"{sys.argv[2]}: order, arcs, R, rank and cover size identical to this run: {same}",flush=True)
okk = Sn==909 and rk==9 and bad==0 and len(nzr)==909 and len(nzc)==909 and cover is not None and len(cover)==906 and s.StatusName(st_)=='OPTIMAL' and not sub.any() and (len(sys.argv)<=2 or same)
print('PASS_EXPXX_SIM909' if okk else 'FAIL_EXPXX_SIM909')
