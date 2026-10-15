"""Dual of the top analysis at dim 1: node <p_s> must reach all 52 tops t with p_s in H_t.
Source s's own dim-2 frame is <p_s,u>, u in K_s. Max union of tops over (u in K_s, y any)."""
import json
from collections import Counter
P=json.load(open('instA.json'))['ports']; cen=json.load(open('cen.json'))['cen']; n=112
nor=[1,7,25,31,121,127,97,103]; dot=lambda a,b:(a&b).bit_count()&1
T=lambda Q: frozenset(i for i,p in enumerate(P) if all(dot(p,q)==0 for q in Q))
best=Counter(); own=Counter(); ex=[]
for s in range(n):
    K=[x for x in range(1,512) if x!=P[s] and all(dot(x,nor[q])==0 for q in cen[s])]
    full=T([P[s]]); assert len(full)==52
    ys=[T([P[s],y]) for y in range(1,512) if y!=P[s]]
    bu=0
    for u in K:
        A=T([P[s],u]); own[len(A)]+=1
        m=max(len(A|B) for B in ys); bu=max(bu,m)
    best[bu]+=1
print('own-L2 top counts over u in K_s (all s):',dict(own))
print('max |reach| with source + ONE other register at <p_s>:',dict(best))
