import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""rho(3): quads X (3-dim, 4 ports = all its odd vectors): r = ports(X^perp) (dim 6) and best k-hyperplane covers of X^perp."""
import json,itertools,sys
from collections import Counter
sys.path.insert(0,'../expBV'); from wordmodel import basis
P=json.load(open('instA.json'))['ports']; n=112; dot=lambda a,b:(a&b).bit_count()&1; isp={p:i for i,p in enumerate(P)}
Mi=[sum(1<<i for i,p in enumerate(P) if dot(p,q)==0) for q in range(512)]
pc=lambda m: bin(m).count('1')
quads=set()
for a,b,c in itertools.combinations(range(n),3):
    d=P[a]^P[b]^P[c]
    if d in isp: quads.add(frozenset((a,b,c,isp[d])))
print('#quads',len(quads))
st=Counter(); best=0
for Q in quads:
    q=list(Q); m=Mi[P[q[0]]]&Mi[P[q[1]]]&Mi[P[q[2]]]; r=pc(m)
    # hyperplanes of X^perp = X^perp cap y^perp for y not in X: ports = m & Mi[y]
    hs=sorted({m&Mi[y] for y in range(1,512)}-{m},key=pc,reverse=True)
    g1=pc(hs[0]); g2=max(pc(x|y) for x in hs[:40] for y in hs[:40])
    ratio=max(1*g1, 4*g2/2, 4*r/3)  # k=1: one pair ->c<=2? (c(1)=2 via one dim2 hyperplane)
    ratio=max(2*g1, 4*g2/2, 4*min(r,3*g1)/3)
    st[(r,g1,g2)]+=1; best=max(best,ratio)
print('rho(3) upper (quads)',best); print(sorted(st.items())[-8:])
