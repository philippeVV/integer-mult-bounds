import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""Price helper inventories (same base and ceiling call as price2.py). argv: JSON list of [label, {width:count}]."""
import sys, math, json; sys.path.insert(0,'../expBV'); sys.path.insert(0,__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', '..', '..', 'complex-local-word-ceiling'))
from wordmodel import Word
from collections import Counter
from fractions import Fraction as Q
import ceiling as CE
W=Word('../../helper-bounds/sel1094.json.gz'); v=W.v; g=Q(1,10**9)
def cap(prof): return CE.ceiling(dict(+Counter(prof)),73,896,9,64,grid=g)[1]
b0=cap(W.paid); p=Counter(W.paid)
for r in W.helpers+list(range(2*v)):
    for c in W.climbs.get(r,[]): p[c[1]]-=1
p[1]+=8*2*v; base=+p
for lab,inv in json.loads(sys.argv[1]):
    q=Counter(base)
    for w,k in inv.items(): q[int(w)]+=k
    print(f'{lab:40s} x{float(cap(q)/b0):.4f} logmass {sum(k*int(w)*math.log(73/int(w)) for w,k in inv.items()):.0f}',flush=True)
