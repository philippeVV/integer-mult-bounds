"""Independent stdlib structural check of schedule{q}.json: arcs are the PG(2,q) Singer incidences (each once), Eulerian k/k,
Gamma - R acyclic by Kahn's algorithm (no use of the stored order), |R| = S(q)."""
import json, sys, math
J=json.load(open(sys.argv[1])); q=J['q']; D=J['D']; N=q*q+q+1; k=(q+1)//2; r=math.isqrt(q)
S=((q+1)*(q*q+q+1)-r*(q*q+1))//4
assert len(set((x-y)%N for x in D for y in D if x!=y))==N-1, "not a perfect difference set"
inc=sorted((i,N+(i+d)%N) for i in range(N) for d in D)
assert sorted(tuple(sorted(a)) for a in J['arcs'])==inc
indeg=[0]*(2*N); outdeg=[0]*(2*N)
for a,b in J['arcs']: outdeg[a]+=1; indeg[b]+=1
assert set(indeg)=={k} and set(outdeg)=={k}
Rs=set(J['R']); adj=[[] for _ in range(2*N)]; deg=[0]*(2*N)
for e,(a,b) in enumerate(J['arcs']):
    if e not in Rs: adj[a].append(b); deg[b]+=1
st=[v for v in range(2*N) if deg[v]==0]; seen=0
while st:
    v=st.pop(); seen+=1
    for w in adj[v]:
        deg[w]-=1
        if deg[w]==0: st.append(w)
print(f"q={q}: incidences OK, Eulerian {k}/{k} OK, Gamma-R acyclic: {seen==2*N}, |R|={len(Rs)}, S(q)={S}, exact count: {len(Rs)==S}")
assert seen==2*N and len(Rs)==S
print("PASS_EXPXQ_JSON")
