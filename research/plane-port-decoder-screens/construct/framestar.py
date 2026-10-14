"""Exact decoder: 8 frame terms (2G) + row-star corrections of every broken odd entry. Verifies sum == 4I on all
non-free pairs (Fractions/ints), prices each term by F2 form rank on plane spans."""
import sys, json; sys.path.insert(0,'.')
from e8 import *
L,labs,G,fr,ports=p112(); n=len(ports)
R=json.load(open(sys.argv[1]))['R']; U=[labs[t] for t in ports]; E=[(U[i],R[i]) for i in range(n)]
def orthp(i,j): return all(gform(a,b,G)==0 for a in E[i] for b in E[j])
FREE=[[i!=j and orthp(i,j) for j in range(n)] for i in range(n)]
Cf=[[ip(L[ports[i]],L[f]) for f in fr] for i in range(n)]
Gm=[[ip(L[ports[i]],L[ports[j]]) for j in range(n)] for i in range(n)]
terms=[]
for k in range(8):
    v={i:Cf[i][k] for i in range(n) if Cf[i][k]}; terms.append((v,v))
for s in range(n):
    B={j:-2*Gm[s][j] for j in range(n) if j!=s and Gm[s][j]%2 and not FREE[s][j]}
    if B: terms.append(({s:1},B))
M=[[0]*n for _ in range(n)]
for a,b in terms:
    for i,x in a.items():
        for j,y in b.items(): M[i][j]+=x*y
bad=sum(1 for i in range(n) for j in range(n) if not FREE[i][j] and M[i][j]!=(4 if i==j else 0))
costs=[frank([w for i in a for w in E[i]],[w for i in b for w in E[i]],G) for a,b in terms]
from collections import Counter
print(sys.argv[1],'bad',bad,'frame costs',costs[:8],'star cost hist',Counter(costs[8:]),'C',sum(costs))
json.dump({'C':sum(costs),'costs':costs,'bad':bad},open(sys.argv[1].replace('.json','_framestar.json'),'w'))
