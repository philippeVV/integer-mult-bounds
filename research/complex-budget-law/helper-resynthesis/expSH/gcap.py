import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""g(a) = max over a-dim A, y notin A of ports(A)*#ports orthogonal to <A,y>  (edge pair-capacity at levels a and 8-a)."""
import json, numpy as np, itertools, sys
sys.path.insert(0,'../expBV'); from wordmodel import basis
P=json.load(open('instA.json'))['ports']; n=112; dot=lambda a,b:(a&b).bit_count()&1
Mi=[sum(1<<i for i,p in enumerate(P) if dot(p,q)==0) for q in range(512)]
def span(B):
    s=[0]
    for v in B: s=s+[x^v for x in s]
    return s
isp=set(P)
from collections import Counter
for a in (1,2,3,4):
    seen=set(); best=Counter()
    # enumerate all a-dim subspaces A spanned by ports (others have fewer ports per dim; ports(A)=0 gives 0)
    # plus subspaces with ports(A)>=1 not port-spanned are dominated? no: include A = span(ports)+free vectors handled by smaller a with free y.
    for S in itertools.combinations(range(n),a):
        B=basis([P[i] for i in S])
        if len(B)<a or B in seen: continue
        seen.add(B); sp=span(B); pa=sum(x in isp for x in sp)
        m=(1<<n)-1
        for b in B: m&=Mi[b]
        spset=set(sp)
        mx=max(bin(m&Mi[y]).count('1') for y in range(1,512) if y not in spset)
        best[pa*mx]+=1
    print('a',a,'#port-spanned A',len(seen),'max cap',max(best),'top values',sorted(best.items())[-5:],flush=True)
