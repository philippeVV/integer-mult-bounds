import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
import json,sys
from collections import Counter,defaultdict
sys.path.insert(0,'../expBV'); from wordmodel import basis
P=json.load(open('instA.json'))['ports']; cen=json.load(open('cen.json'))['cen']; n=112
nor=[1,7,25,31,121,127,97,103]
dot=lambda a,b:(a&b).bit_count()&1
cl=defaultdict(list)
for s in range(n): cl[tuple(cen[s])].append(s)
C=list(cl.values()); cid={s:i for i,c in enumerate(C) for s in c}
def span(B):
    s=[0]
    for v in B: s=s+[x^v for x in s]
    return s
X=[basis([P[i] for i in c]) for c in C]
anti=[[j for j in range(14) if sum(dot(P[a],P[b])==0 for a in C[i] for b in C[j])==0] for i in range(14)]
print('antipodes',anti)
# intersections of class spaces
print('dim X_i cap X_j', Counter(len(basis(X[i]))+len(basis(X[j]))-len(basis(list(X[i])+list(X[j]))) for i in range(14) for j in range(i+1,14)))
# within class: orth graph structure
i=0; c=C[0]
print('own-class orth graph degrees',[sum(dot(P[a],P[b])==0 for b in c if b!=a) for a in c])
print('own-class orth: is complement of pairing?', [[c.index(b) for b in c if b!=a and dot(P[a],P[b])==0] for a in c])
# between class 0 and class j: biadjacency structure
for j in range(1,14):
    if j in anti[0]: continue
    M=[[int(dot(P[a],P[b])==0) for b in C[j]] for a in C[0]]
    rows=Counter(tuple(r) for r in M)
    print(j,'distinct rows',len(rows),'row mult',sorted(rows.values()), 'sum of orth a: ', Counter(P[a]^P[b] for a in C[0] for b in C[j] if dot(P[a],P[b])==0).most_common(2)[:1])
    if j==1: print(M)
# normals vs classes
print('normal_q in class spaces?',[[int(len(basis(list(X[i])+[nor[q]]))==4) for q in range(8)] for i in range(3)])
print('class sizes',[len(c) for c in C], 'cen of classes',list(cl.keys()))
json.dump(dict(C=C,anti=anti),open('classes.json','w'))
