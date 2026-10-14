import itertools
from fractions import Fraction as F
# E8 roots: D8 roots (+-e_i +-e_j) and (+-1/2)^8 with even number of minus signs
R=[]
for i,j in itertools.combinations(range(8),2):
    for si in (1,-1):
        for sj in (1,-1):
            v=[0]*8; v[i]=si; v[j]=sj; R.append(tuple(F(x) for x in v))
for signs in itertools.product((1,-1),repeat=8):
    if signs.count(-1)%2==0: R.append(tuple(F(s,2) for s in signs))
assert len(R)==240
ip=lambda a,b: sum(x*y for x,y in zip(a,b))
lines=[]; seen=set()
for r in R:
    k=min(r,tuple(-x for x in r))
    if k not in seen: seen.add(k); lines.append(k)
assert len(lines)==120
# a frame: 8 mutually orthogonal lines: (e0+-e1),(e2+-e3),... 
frame=[i for i,l in enumerate(lines) if sum(1 for x in l if x!=0)==2 and [k for k,x in enumerate(l) if x!=0] in ([0,1],[2,3],[4,5],[6,7])]
assert len(frame)==8 and all(ip(lines[a],lines[b])==0 for a,b in itertools.combinations(frame,2))
P=[i for i in range(120) if i not in frame]; n=len(P)
odd=[[abs(ip(lines[P[a]],lines[P[b]]))==1 for b in range(n)] for a in range(n)]
print('P112 odd degree', set(sum(r) for r in odd))
nb=[set(b for b in range(n) if odd[a][b]) for a in range(n)]
K5=[]
def ext(cl,cand):
    if len(cl)==5: K5.append(cl); return
    for v in sorted(cand):
        if not cl or v>cl[-1]: ext(cl+[v],cand&nb[v])
ext([],set(range(n)))
print('K5 count',len(K5))
from pysat.solvers import Cadical153
s=Cadical153()
for k in K5:
    s.add_clause([a+1 for a in k]); s.add_clause([-(a+1) for a in k])
print('2-colouring without monochromatic odd K5 exists:', s.solve())
