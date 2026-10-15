import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
import json, sys
from collections import Counter, defaultdict
sys.path.insert(0,'../expBV'); from wordmodel import basis
I=json.load(open('instA.json')); P=I['ports']; C=I['centres']
U=[tuple(basis(c['U'])) for c in C]
nor=[1,7,25,31,121,127,97,103]
dot=lambda a,b:(a&b).bit_count()&1
cen=[[q for q in range(8) if dot(nor[q],P[s])==0] for s in range(112)]
assert all(sorted(cen[s])==sorted(q for q in range(8) if str(s) in C[q]['src'] or s in map(int,C[q]['src'])) for s in range(112))
Nsp=basis(nor); print('dim span normals',len(Nsp))
# syndrome = pattern of dots with the 8 normals
cls=Counter(tuple(cen[s]) for s in range(112)); print('#classes',len(cls),'sizes',sorted(cls.values()))
def meet_dim(qs):
    return 9-len(basis([nor[q] for q in qs]))
print('K_s dims',Counter(meet_dim(cen[s]) for s in range(112)))
K=[tuple(basis([x for x in range(1,512) if all(dot(x,nor[q])==0 for q in cen[s])])) for s in range(112)]
inK=lambda p,k: basis(list(K[k])+[p])==K[k]
rel=Counter((tuple(cen[a])==tuple(cen[b]), inK(P[a],b)) for a in range(112) for b in range(112) if a!=b)
print('(same class, p_a in K_b):',rel)
# distinct port values?
print('distinct ports',len(set(P)))
# class spans
for k,c in list(cls.items())[:3]:
    ss=[s for s in range(112) if tuple(cen[s])==k]; print(k,len(ss),'span dim',len(basis([P[s] for s in ss])))
json.dump(dict(cen=cen),open('cen.json','w'))
