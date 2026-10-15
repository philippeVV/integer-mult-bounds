import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
import sys, math; sys.path.insert(0,'../expBV')
from wordmodel import Word, basis, regs
from fractions import Fraction as F
from collections import Counter, defaultdict
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
for sec in ('precondition','compensation','A','cut'): 
    for e in w[sec]: ap(e)
aftercut=[dict(r) for r in rows]
for e in w['B']: ap(e)
afterB=[dict(r) for r in rows]
for e in w['restore']: ap(e)
# target rows (source columns only, helper columns = dirt)
ok=0; sz=Counter(); spandim=Counter(); orth=0
for t in range(v):
    r=afterB[v+t]; src={k:x for k,x in r.items() if k<v}
    if src=={t:F(1)}: ok+=1
    c={k:x for k,x in aftercut[v+t].items() if k<v}
    corr={k: (1 if k==t else 0)-c.get(k,0) for k in set(c)|{t}}
    corr={k:x for k,x in corr.items() if x}
    sz[len(corr)]+=1
    spandim[len(basis([ports[k] for k in corr]))]+=1
    pt=ports[t]
    orth+= all((ports[k]&pt).bit_count()%2==0 for k in corr)
    if t<3: print(t, 'cut-row src',len(c), 'corr',len(corr), Counter(corr.values()))
print('afterB target src rows exact',ok,'/',v)
print('corr support sizes',sorted(sz.items()),'span dims',sorted(spandim.items()),'all corr ports ⊥ p_t',orth)
# p.p parity
print('port self-parity', Counter((p).bit_count()%2 for p in ports))
# sources rows after B / restore
print('sources restored after restore:', sum(1 for s in range(v) if {k:x for k,x in rows[s].items() if k<2*v}=={s:F(1)}))
print('targets exact after restore:', sum(1 for t in range(v) if {k:x for k,x in rows[v+t].items() if k<2*v}=={t:F(1)}))
