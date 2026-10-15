import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""Price scenarios for the selected #1094 word with the repo ceiling tool (same call as propdag.py)."""
import sys, math, json; sys.path.insert(0,'../expBV'); sys.path.insert(0,__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', '..', '..', 'complex-local-word-ceiling'))
from wordmodel import Word
from collections import Counter
from fractions import Fraction as Q
import ceiling as CE
W=Word('../../helper-bounds/sel1094.json.gz'); v=W.v
g=Q(1,10**12)
def cap(prof): return CE.ceiling(dict(+Counter(prof)),73,896,9,64,grid=g)[1]
b0=cap(W.paid)
def strip(regs):
    p=Counter(W.paid)
    for r in regs:
        for c in W.climbs.get(r,[]): p[c[1]]-=1
    return p
data=list(range(2*v))
out={}
p=strip(W.helpers); out['ideal_no_helper_climbs']=p
q=Counter(p); q=strip(W.helpers+data); q[1]+=8*len(data); out['no_helpers_data_all_rank1']=q
q=strip(W.helpers+data); q[8]+=len(data); out['no_helpers_data_oneshot']=q
# helpers one-shot at their pinned widths (any gate multiset with this inventory)
q=strip(W.helpers)
for r in W.helpers:
    w=len(W.end[r])-len(W.start[r])
    if w: q[w]+=1
out['helpers_oneshot_fixed_inventory']=q
# minimal helper inventory from first-touch lemma: 112 helpers with a rank-1 climb only, data worst
q=strip(W.helpers+data); q[1]+=8*len(data)+112; out['firsttouch_112_rank1_data_worst']=q
res={}
for k,pr in out.items():
    b=cap(pr); res[k]=float(b/b0)
    print(f'{k:40s} ceiling x{float(b/b0):.4f}  logmass {sum(n*r*math.log(73/r) for r,n in pr.items()):.0f}')
print('baseline logmass',round(sum(n*r*math.log(73/r) for r,n in W.paid.items())))
json.dump(res,open('price.json','w'),indent=1)
