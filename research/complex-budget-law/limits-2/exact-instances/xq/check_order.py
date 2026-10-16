"""Independent stdlib checker: usage check_order.py graph.txt order.txt [S-split t]
Graph file: 'V d' then V lines of neighbours. Prints Boot cost sum_v (d/2 - j_v)^+, max prefix cut,
total ascent/2 of the cut profile (must equal the cost), and the schedule back-arc flow feasibility (boot3 rule)."""
import sys
from collections import deque
L=open(sys.argv[1]).read().split('\n'); h=L[0].split(); V=int(h[0]); d=int(h[1]) if len(h)>1 else 16; K=d//2
nb=[list(map(int,L[1+v].split())) for v in range(V)]
assert all(len(x)==d for x in nb) and all(v in nb[u] for v in range(V) for u in nb[v])
order=[int(x) for x in open(sys.argv[2]).read().split()]; assert sorted(order)==list(range(V))
pos=[0]*V
for i,v in enumerate(order): pos[v]=i
j=[sum(1 for u in nb[v] if pos[u]<pos[v]) for v in range(V)]
cost=sum(max(0,K-x) for x in j)
c=0; prof=[]; asc=0
for v in order:
    dc=d-2*j[v]; asc+=max(0,dc); c+=dc; prof.append(c)
assert c==0 and asc==2*cost
print('V',V,'deg',d,'Boot cost',cost,'max prefix cut',max(prof),'peaks-check ok')
# back-arc flow: excess vertices (j>K) send j-K back arcs to earlier deficit nbrs; Dinic-free simple augmenting (unit caps)
if len(sys.argv)>3 and sys.argv[3]=="noflow": sys.exit(0)
need=sum(max(0,x-K) for x in j)
src=[v for v in range(V) if j[v]>K]
cap={}
G={}
def add(a,b,cp):
    G.setdefault(a,[]).append(b); G.setdefault(b,[]).append(a); cap[(a,b)]=cap.get((a,b),0)+cp; cap.setdefault((b,a),0)
for v in src:
    add('s',('o',v),j[v]-K)
    for u in nb[v]:
        if pos[u]<pos[v] and j[u]<K: add(('o',v),('i',u),1)
for u in range(V):
    if j[u]<K: add(('i',u),'t',K-j[u])
flow=0
while True:
    par={'s':None}; dq=deque(['s'])
    while dq and 't' not in par:
        a=dq.popleft()
        for b in G.get(a,[]):
            if b not in par and cap[(a,b)]>0: par[b]=a; dq.append(b)
    if 't' not in par: break
    b='t'; f=10**9
    while par[b] is not None: f=min(f,cap[(par[b],b)]); b=par[b]
    b='t'
    while par[b] is not None: cap[(par[b],b)]-=f; cap[(b,par[b])]+=f; b=par[b]
    flow+=f
print('back-arc flow',flow,'/',need,'schedule from this order:', flow==need)
