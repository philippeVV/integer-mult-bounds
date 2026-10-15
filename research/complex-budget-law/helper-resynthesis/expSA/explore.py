import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
import sys; sys.path.insert(0,'../expBV')
from wordmodel import Word, basis, regs
from fractions import Fraction as F
from collections import Counter, defaultdict
W=Word('../../helper-bounds/sel1094.json.gz'); w=W.w; v=W.v; n=W.n
def role(r): return 'S' if r<v else ('T' if r<2*v else 'H')
print('centers',w['centers'], [len(f) for f in w['center_frames']])
# section x role-pair counts and frame dims
for sec in ('precondition','compensation','A','cut','B','restore'):
    c=Counter()
    for e in w[sec]:
        if e[0]!='add': c[('scale',role(e[1]))]+=1; continue
        U=e[-1]; d=None if U is None else len(basis(U))
        c[(role(e[1])+'+='+role(e[2]),d)]+=1
    print(sec,len(w[sec]),sorted(c.items(),key=lambda x:-x[1])[:14])
# scalar rows after precondition+compensation+A
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
ports=w['ports']
for q,U in zip(w['centers'],w['center_frames']):
    r=rows[q]; src={k:x for k,x in r.items() if k<v}; hel={k:x for k,x in r.items() if k>=2*v}
    print('centre',q,'port',ports[q] if q<v else None,'U',U,'#src',len(src),'coef',Counter(src.values()),'#hel',len(hel))
    Ub=basis(U)
    ok=all(basis(list(Ub)+[ports[k]])==Ub for k in src)
    print('   src ports in U:',ok)
# cut structure
cc=Counter((e[1]-v,e[2]) for e in w['cut'])
pert=Counter(t for t,q in cc)
print('cut per target',Counter(pert.values()), 'targets',len(pert))
