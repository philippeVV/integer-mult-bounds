import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
import json, sys
from collections import Counter
sys.path.insert(0,'../expBV'); from wordmodel import basis
P=json.load(open('instA.json'))['ports']; n=112; dot=lambda a,b:(a&b).bit_count()&1
idx={p:i for i,p in enumerate(P)}
S=lambda Q: frozenset(i for i,p in enumerate(P) if all(dot(p,q)==0 for q in Q))
G={t:[u for u in range(n) if u!=t and len(S([P[t],P[u]]))==26] for t in range(n)}
print('G26 degree',Counter(len(G[t]) for t in range(n)))
print('t,u in G26: p_t.p_u parity',Counter(dot(P[t],P[u]) for t in range(n) for u in G[t]), ' p_t+p_u is port',Counter((P[t]^P[u]) in idx for t in range(n) for u in G[t]))
pairs={t:[(a,b) for i,a in enumerate(G[t]) for b in G[t][i+1:] if len(S([P[t],P[a]])|S([P[t],P[b]]))==52] for t in range(n)}
print('valid neighbour pairs per t',Counter(len(pairs[t]) for t in range(n)))
t=0; print('example t=0 G',G[0][:20],'pairs',pairs[0][:10])
# structure: is the valid-pair relation 'b = a + something'?
print('p_a+p_b+p_t port?',Counter((P[a]^P[b]^P[t]) in idx for t in range(n) for a,b in pairs[t]), ' p_a.p_b',Counter(dot(P[a],P[b]) for t in range(n) for a,b in pairs[t]))
json.dump(dict(G=G,pairs={t:pairs[t] for t in pairs}),open('g26.json','w'))
