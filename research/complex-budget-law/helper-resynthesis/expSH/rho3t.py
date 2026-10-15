import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""rho(3) = rho(6) <= 24. Node X (dim 3): ratio = max_k f(k)g(k)/k.  k=1 term <= E25 = 24 (e25.py); k>=3 term <= 4*20/3 < 26.
k=2 term = f(2)*g(2)/2 with f(2)<=ports(X)<=4: quads checked in rho3.py (g2<=12 -> 24). Here: triples (3 ports): 3*g2/2 <= 3*r/2,
so only triples with r=ports(X^perp)>=17 need g2; non-port-spanned X have in=1 (k=1 term, <= E25)."""
import json,itertools,sys
sys.path.insert(0,'../expBV'); from wordmodel import basis
P=json.load(open('instA.json'))['ports']; n=112; dot=lambda a,b:(a&b).bit_count()&1; isp={p:i for i,p in enumerate(P)}
Mi=[sum(1<<i for i,p in enumerate(P) if dot(p,q)==0) for q in range(512)]
pc=lambda m: bin(m).count('1')
from collections import Counter
st=Counter(); worst=0
for a,b,c in itertools.combinations(range(n),3):
    if P[a]^P[b]^P[c] in isp: continue
    m=Mi[P[a]]&Mi[P[b]]&Mi[P[c]]; r=pc(m); st[r]+=1
    if r>=17:
        sp={0,P[a],P[b],P[c],P[a]^P[b],P[a]^P[c],P[b]^P[c],P[a]^P[b]^P[c]}
        hs=list({m&Mi[y] for y in range(1,512) if y not in sp}); g2=max(pc(x|y) for x in hs for y in hs); worst=max(worst,3*g2/2)
print('triples: ports(X^perp) distribution',sorted(st.items()),' worst k=2 ratio among r>=17:',worst)
print('rho(3)=rho(6) <=',max(24,worst,3*16/2))
