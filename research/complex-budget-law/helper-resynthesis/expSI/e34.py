import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""E34 = max ports(A)*ports(B), A perp B, dim A<=3, dim B<=4. Enumerate port-spanned 4-dim B (all), max ports in 3-dim subspace of B^perp."""
import json, numpy as np, itertools, sys
sys.path.insert(0,'../expBV'); from wordmodel import basis
P=json.load(open('instA.json'))['ports']; n=112; dot=lambda a,b:(a&b).bit_count()&1; isp={p:i for i,p in enumerate(P)}
Mi=[sum(1<<i for i,p in enumerate(P) if dot(p,q)==0) for q in range(512)]
def span(B):
    s=[0]
    for v in B: s=s+[x^v for x in s]
    return s
def maxports3(idx):
    # max #ports of the list lying in a 3-dim subspace: quads (sum zero 4-sets) -> 4, else min(3,len)
    if len(idx)<=3: return len(idx)
    S=set(P[i] for i in idx)
    for a,b,c in itertools.combinations(idx,3):
        if P[a]^P[b]^P[c] in S: return 4
    return 3
seen=set(); best=0; from collections import Counter; st=Counter()
for S4 in itertools.combinations(range(n),4):
    B=basis([P[i] for i in S4])
    if len(B)<4 or B in seen: continue
    seen.add(B)
    pb=sum(1 for x in span(B) if x in isp)
    m=(1<<n)-1
    for b in B: m&=Mi[b]
    perp=[i for i in range(n) if m>>i&1]
    v=pb*maxports3(perp); st[(pb,len(perp))]+=1
    if v>best: best=v; print('new',v,pb,perp,flush=True)
print('E34',best,'#B',len(seen)); print(sorted(st.items()))
