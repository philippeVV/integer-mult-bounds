# negative controls for check_F0 case B encoding
import check_F0 as C
from pysat.solvers import Cadical153
import itertools
x=2
# (1) drop FR clauses -> ?
src=open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'check_F0.py')).read()
def run(skip_fr=False, skip_tri=False, only_first_tri=False):
    Sx=set(r for r in C.LIVE if C.ip(x,r)==0 or r==x); nv=[0]
    def new(): nv[0]+=1; return nv[0]
    a={r:new() for r in C.LIVE}; b={r:new() for r in C.LIVE}; S=Cadical153(); add=S.add_clause
    for r in C.LIVE:
        if r not in Sx: add([-a[r]]); add([-b[r]])
    if not skip_fr:
        for i in range(28):
            for j in range(28):
                if i!=j and not (C.CL[i][0]&C.CL[j][0]):
                    for r in C.CL[i][1]:
                        for s in C.CL[j][1]: add([-a[r],-b[s]])
    pl=[]
    for sv,Q in C.CL:
        na,nb,p=new(),new(),new()
        for r in Q: add([-na,-a[r]]); add([-nb,-b[r]])
        add([na]+[a[r] for r in Q]); add([nb]+[b[r] for r in Q]); add([-p,na,nb]); add([p,-na]); add([p,-nb]); pl.append(p)
    ab={}
    def AB(r,s):
        if (r,s) not in ab:
            t=new(); add([-t,a[r]]); add([-t,b[s]]); add([t,-a[r],-b[s]]); ab[r,s]=t
        return ab[r,s]
    if not skip_tri:
        for P,*parts in C.TRI:
            xs=[]
            for (i,X,Xb),(j,Y,Yb) in itertools.combinations(parts,2):
                for U,W in ((X,Y),(Y,X),(Xb,Yb),(Yb,Xb)): xs+=[AB(r,s) for r in U for s in W]
            ids=[p[0] for p in parts]
            for k in range(3):
                o=[ids[m] for m in range(3) if m!=k]; add([-pl[o[0]],-pl[o[1]]]+xs)
    r=S.solve()
    if r:
        m=set(l for l in S.get_model() if l>0)
        return r, sum(a[q] in m for q in C.LIVE), sum(b[q] in m for q in C.LIVE)
    return r
print('full', run())
print('no FR', run(skip_fr=True))
print('no TRI', run(skip_tri=True))
