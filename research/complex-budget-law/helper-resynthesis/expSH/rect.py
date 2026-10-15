import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""Pair-rectangle bound: a delivered pair (s,t) passes, at every dim d, a level L with p_s in L and p_t in L^perp.
Max ports(L)*ports(L^perp) over dim-d subspaces spanned (WLOG) by ports. Random+greedy search per d."""
import json, random, sys
sys.path.insert(0,'../expBV'); from wordmodel import basis
P=json.load(open('instA.json'))['ports']; dot=lambda a,b:(a&b).bit_count()&1
def span(vs):
    s={0}
    for v in vs: s|={x^v for x in s}
    return s
def score(B):
    sp=span(B); a=sum(p in sp for p in P); b=sum(all(dot(p,x)==0 for x in B) for p in P); return a,b
random.seed(0); best={}
for d in range(1,9):
    bb=(0,0,0)
    for trial in range(3000):
        B=[]; 
        while len(B)<d:
            # grow by ports orthogonal-friendly choice
            c=random.choice(P) if random.random()<0.85 else random.randrange(1,512)
            if len(basis(B+[c]))==len(B)+1: B.append(c)
        # local improvement: replace one vector
        a,b=score(B); cur=a*b
        for _ in range(60):
            i=random.randrange(d); c=random.choice(P) if random.random()<0.8 else random.randrange(1,512)
            B2=B[:i]+[c]+B[i+1:]
            if len(basis(B2))<d: continue
            a2,b2=score(B2)
            if a2*b2>=cur: B,cur,a,b=B2,a2*b2,a2,b2
        if cur>bb[0]: bb=(cur,a,b)
    best[d]=bb; print('d',d,'max ports(L)*ports(L^perp) found',bb, ' => levels needed >=',-(-5824//bb[0]) if bb[0] else None, flush=True)
