import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
import json, sys
from collections import Counter
sys.path.insert(0,'../expBV'); from wordmodel import basis
P=json.load(open('instA.json'))['ports']; n=112; dot=lambda a,b:(a&b).bit_count()&1
idx={p:i for i,p in enumerate(P)}
g=json.load(open('g26.json')); G={int(k):v for k,v in g['G'].items()}
S=lambda Q: frozenset(i for i,p in enumerate(P) if all(dot(p,q)==0 for q in Q))
def span(vs):
    s={0}
    for v in vs: s|={x^v for x in s}
    return s
stats=Counter(); cov=Counter(); sh=Counter()
for (t,a) in [(0,G[0][0]),(0,G[0][1]),(5,G[5][3])]:
    D=span([P[t],P[a]]); full=S([P[t],P[a]])
    hyp={}
    for c in range(1,512):
        if c in D: continue
        key=tuple(sorted(span([P[t],P[a],c])))
        if key in hyp: continue
        hyp[key]=(c,S([P[t],P[a],c]))
    sizes=Counter(len(v[1]) for v in hyp.values()); print('edge',t,a,'#hyperplanes',len(hyp),'port counts',sorted(sizes.items()))
    hl=list(hyp.items())
    covers=[(k1,k2) for i,(k1,v1) in enumerate(hl) for k2,v2 in hl[i+1:] if len(v1[1]|v2[1])==26]
    print('  2-covers of the 26 ports:',len(covers), Counter((len(hyp[k1][1]),len(hyp[k2][1])) for k1,k2 in covers))
    # ports in each dual plane (other tops sharing L6)
    for k1,k2 in covers[:3]:
        for k in (k1,k2):
            pp=[idx[x] for x in k if x in idx]; print('   L6 dual plane ports',pp,'G26 edges inside',[(x,y) for x in pp for y in pp if x<y and y in G[x]])
