import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
import sys, math; sys.path.insert(0,'../expBV')
from wordmodel import Word, basis, regs
from collections import Counter, defaultdict
W=Word('../../helper-bounds/sel1094.json.gz'); w=W.w; v=W.v
f=lambda r: r*math.log(73/r)
secof=lambda t: W.events[t][0] if t<W.T_END else 'END'
for kind,rr in (('helper',W.helpers),('source',range(v)),('target',range(v,2*v))):
    m=Counter(); lm=Counter(); 
    for r in rr:
        for (t,d,a,b) in W.climbs.get(r,[]):
            m[secof(t)]+=d; lm[secof(t)]+=f(d)
    print(kind, 'mass',dict(m),'log',{k:round(x) for k,x in lm.items()}, 'tot',round(sum(lm.values())))
# helper touch sections pattern
pat=Counter()
for r in W.helpers:
    secs=[]
    for t,s in W.touches[r]:
        s=s.split(':')[0]
        if s in('precondition','compensation'): s='pre'
        if not secs or secs[-1]!=s: secs.append(s)
    pat[tuple(secs)]+=1
for k,c in pat.most_common(15): print(c,k)
