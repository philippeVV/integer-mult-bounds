"""Is the per-port charge-3 bound attainable locally? For each s: N_s = same-class pair node (charge 1) + one 1-port node <p_s,y>
(charge 2, reach<=32) with reach union = 52; or N_s 1-port (u in K_s non-port) + G26 pair etc."""
import json
from collections import Counter
P=json.load(open('instA.json'))['ports']; cen=json.load(open('cen.json'))['cen']; n=112; isp={p:i for i,p in enumerate(P)}
dot=lambda a,b:(a&b).bit_count()&1
T=lambda Q: frozenset(i for i,p in enumerate(P) if all(dot(p,q)==0 for q in Q))
res=Counter()
for s in range(n):
    full=T([P[s]])
    one=[T([P[s],y]) for y in range(1,512) if y!=P[s] and y not in isp and (y^P[s]) not in isp]
    g26=[T([P[s],P[q]]) for q in range(n) if q!=s and len(T([P[s],P[q]]))==26]
    a=any(len(T([P[s],P[q]])|B)==52 for q in range(n) if q!=s and cen[q]==cen[s] for B in one)
    Ks=[u for u in range(1,512) if u!=P[s] and u not in isp and (u^P[s]) not in isp and all(dot(u,[1,7,25,31,121,127,97,103][k])==0 for k in cen[s])]
    b=any(len(T([P[s],u])|B|C)==52 for u in Ks for i,B in enumerate(g26) for C in g26[i+1:])
    res[(a,b)]+=1
print('(same-class N_s + 1-port node covers 52, 1-port N_s + two G26 covers 52):',dict(res))
