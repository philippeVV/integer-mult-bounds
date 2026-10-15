import os
import pickle,sys; sys.path.insert(0,'.')
from frames import *
from price_feas import price,W,feas
best={}
for c in feas: best.setdefault(c['h'],c)
for hs in ([606,669],[606,638,639,669,670,671],list(best)):
    ns={h:best[h]['S'] for h in hs}
    fr,fl,t=W.lift(ns); sd=Counter()
    for r in t:
        for d in W.climbs(r): sd[d]-=1
        for d in W.climbs(r,fr,ns.get(r)): sd[d]+=1
    sd={k:v for k,v in sd.items() if v}
    print(hs,sorted(fl),sd,'%.3e'%price(sd,[(5,4)]*len(hs)))
