import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""Price helper inventories on top of the zero-helper-mass base (data registers charged worst case: 8 rank-1
climbs each), same ceiling call as expSA/price.py. One-shot helpers (each helper = one climb of rank = width)."""
import sys, math, json; sys.path.insert(0,'../expBV'); sys.path.insert(0,__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', '..', '..', 'complex-local-word-ceiling'))
from wordmodel import Word
from collections import Counter
from fractions import Fraction as Q
import ceiling as CE
W=Word('../../helper-bounds/sel1094.json.gz'); v=W.v; g=Q(1,10**9)
def cap(prof): return CE.ceiling(dict(+Counter(prof)),73,896,9,64,grid=g)[1]
b0=cap(W.paid)
p=Counter(W.paid)
for r in W.helpers+list(range(2*v)):
    for c in W.climbs.get(r,[]): p[c[1]]-=1
p[1]+=8*2*v; base=+p
res={}
def run(lab,extra):
    q=Counter(base)
    for w,k in extra.items(): q[w]+=k
    x=float(cap(q)/b0); res[lab]=x; print(f'{lab:34s} x{x:.4f}  helper logmass {sum(k*w*math.log(73/w) for w,k in extra.items()):.0f}',flush=True)
run('zero helpers',{})
run('8 x w1 (centre lower bound)',{1:8})
run('19 x w1 (dim-6 deficit 224/12)',{1:19})
run('112 x w1',{1:112})
run('112 x w2',{2:112})
run('112 x w8',{8:112})
run('240 x w6',{6:240})
run('240 x w9',{9:240})
run('today helper inventory one-shot',{4:56,5:32,7:30,8:122,9:240})
# break-even: number of width-9 / width-1 helpers giving x1.5
for w in (1,3,9):
    lo,hi=0,4000
    while hi-lo>1:
        mid=(lo+hi)//2; q=Counter(base); q[w]+=mid
        if cap(q)/b0>1.5: lo=mid
        else: hi=mid
    print(f'break-even x1.5 with width-{w} helpers: {hi} helpers (logmass {hi*w*math.log(73/w):.0f})',flush=True); res[f'breakeven_w{w}']=hi
json.dump(res,open('price2.json','w'),indent=1)
