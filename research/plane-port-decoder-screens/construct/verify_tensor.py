"""Independent exact check of a tensor(2-class)+private plane-port decoder.
input json: cls, R (r_t bitmask per port index), drop. Builds terms explicitly (integers),
checks: every plane nondegenerate; sum of terms == 4*I on every NON-free pair (free = planes orthogonal);
computes per-term cost = F2 rank of the label form on (sum of row planes) x (sum of col planes)."""
import sys, json; sys.path.insert(0,'.')
from e8 import *
L,labs,G,fr,ports=p112(); n=len(ports)
sol=json.load(open(sys.argv[1])); cls,R,drop=sol['cls'],sol['R'],sol.get('drop',[0]*n)
U=[labs[t] for t in ports]
E=[(U[i],R[i]) for i in range(n)]
for i in range(n):
    gm=[[gform(a,b,G) for b in E[i]] for a in E[i]]
    assert (gm[0][0]*gm[1][1]+gm[0][1]*gm[1][0])%2==1, ('degenerate plane',i)
def orthp(i,j): return all(gform(a,b,G)==0 for a in E[i] for b in E[j])
Cf=[[ip(L[ports[i]],L[f]) for f in fr] for i in range(n)]
terms=[]  # (rowvec dict, colvec dict)
for k in range(8):
    for j in range(2):
        v={i:Cf[i][k] for i in range(n) if not drop[i] and cls[i]==j and Cf[i][k]}
        if v: terms.append((v,dict(v)))
for i in range(n):
    if drop[i]: terms.append(({i:2},{i:2}))
M=[[0]*n for _ in range(n)]
for a,b in terms:
    for i,x in a.items():
        for j,y in b.items(): M[i][j]+=x*y
bad=0; free=0
for i in range(n):
    for j in range(n):
        if i!=j and orthp(i,j): free+=1; continue
        if M[i][j]!=(4 if i==j else 0): bad+=1
costs=[frank([w for i in a for w in E[i]],[w for i in b for w in E[i]],G) for a,b in terms]
print('terms',len(terms),'bad',bad,'free ordered',free,'costs',costs,'C',sum(costs),'dropped',sum(drop))
