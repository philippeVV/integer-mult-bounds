import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""Enumerate all 4-dim X (with an odd vector) : X = E + <o>, E 3-dim even. Record ports(X), ports(X^perp)."""
import json, sys, pickle
from collections import Counter
sys.path.insert(0,'../expBV'); from wordmodel import basis
P=json.load(open('instA.json'))['ports']; n=112; isp={p:i for i,p in enumerate(P)}
nor=[1,7,25,31,121,127,97,103]
dot=lambda a,b:(a&b).bit_count()&1
Mi=[sum(1<<i for i,p in enumerate(P) if dot(p,q)==0) for q in range(512)]
ev=[x for x in range(1,512) if x.bit_count()%2==0]
seen=set(); hubs=[]
for i,a in enumerate(ev):
  for j in range(i+1,len(ev)):
    b=ev[j]; ab=a^b
    for k in range(j+1,len(ev)):
      c=ev[k]
      if c==ab or c==a^b: pass
      if c in (ab,): continue
      if c^a==b: continue
      E=basis([a,b,c])
      if len(E)<3 or E in seen: continue
      seen.add(E)
      el=[0,a,b,ab,c,c^a,c^b,c^ab]
      key={}
      for s,p in enumerate(P):
        kk=min(p^e for e in el); key.setdefault(kk,[]).append(s)
      mE=Mi[a]&Mi[b]&Mi[c]
      for kk,ss in key.items():
        if len(ss)<3: continue
        perp=mE&Mi[kk]
        npp=bin(perp).count('1')
        if npp<3: continue
        hubs.append((E,kk,tuple(ss),perp))
print('3-dim even subspaces',len(seen),'hubs kept',len(hubs))
st=Counter((len(h[2]),bin(h[3]).count('1')) for h in hubs)
for k in sorted(st): print(k,st[k])
pickle.dump(hubs,open('hubs.pkl','wb'))
