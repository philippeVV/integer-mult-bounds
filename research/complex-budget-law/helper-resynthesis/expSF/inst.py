import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""Extract the A-phase centre instance from the selected #1094 word: ports, centre frames, centre supports."""
import sys, json; sys.path.insert(0,'../expBV')
from wordmodel import Word, basis
from fractions import Fraction as F
W=Word('../../helper-bounds/sel1094.json.gz'); w=W.w; v=W.v; n=W.n; ports=w['ports']
rows=[{k:F(1)} for k in range(n)]
def ap(e):
    if e[0]=='add':
        t,s=e[1],e[2]; q=F(e[3],e[4])
        for k,x in list(rows[s].items()):
            a=rows[t].get(k,0)+q*x
            if a: rows[t][k]=a
            else: rows[t].pop(k,None)
    else:
        q=F(e[2],e[3]); rows[e[1]]={k:q*x for k,x in rows[e[1]].items()}
for sec in ('precondition','compensation','A'):
    for e in w[sec]: ap(e)
C=[]
for q,U in zip(w['centers'],w['center_frames']):
    C.append(dict(reg=q,U=list(basis(U)),src={k:int(x) for k,x in rows[q].items() if k<v}))
json.dump(dict(ports=ports,centres=C),open('instA.json','w'))
# normals of centre hyperplanes
def normal(U):
    return [c for c in range(1,512) if all((c&u).bit_count()%2==0 for u in U)]
print([normal(c['U']) for c in C])
from collections import Counter
inc=Counter(k for c in C for k in c['src']); print('sources per #centres',Counter(inc.values()))
