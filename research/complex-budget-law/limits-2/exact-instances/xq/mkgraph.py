"""PG(2,q) Singer graph, q=p^2: point i ~ line N+(i+d)%N, d in D.  writes graph{q}.txt and D{q}.txt. usage: mkgraph.py p"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # python -I omits the script dir
from field import *
open(f'D{q}.txt','w').write(' '.join(map(str,D))+'\n')
nb=[[] for _ in range(2*N)]
for i in range(N):
    for d in D:
        l=(i+d)%N; nb[i].append(N+l); nb[N+l].append(i)
with open(f'graph{q}.txt','w') as f:
    f.write(f'{2*N} {q+1}\n')
    for v in range(2*N): f.write(' '.join(map(str,nb[v]))+'\n')
print(q,N,D)
