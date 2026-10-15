import json, itertools, random
from collections import Counter
P=json.load(open('instA.json'))['ports']; dot=lambda a,b:(a&b).bit_count()&1
# max #ports in Q^perp for Q of dim k (k=1..4), exhaustive for k<=2, sampled for 3,4
def cnt(Q): return sum(all(dot(p,q)==0 for q in Q) for p in P)
print('k=1 max/min',max(cnt([q]) for q in range(1,512)),min(cnt([q]) for q in range(1,512)))
c2=Counter(); m2=0
for a in range(1,512):
    for b in range(a+1,512):
        c=cnt([a,b]); m2=max(m2,c)
print('k=2 max',m2)
# for target t: 52 ports in p_t^perp; max inside a 7-dim subspace = max_a cnt([p_t,a])
print('per target max L7 ports',Counter(max(cnt([P[t],a]) for a in range(1,512) if a!=P[t]) for t in range(112)))
print('per target L7 = <p_t,p_t2>^perp ports',Counter(cnt([P[t],P[u]]) for t in range(112) for u in range(112) if u!=t))
