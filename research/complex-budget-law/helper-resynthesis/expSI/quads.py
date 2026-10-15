import sys as _s, os as _o; _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "..", "helper-bounds"))  # HB_PATH
"""All ordered orthogonal pairs (s,t) (+ centre pairs (s,q), p_s . nor_q = 0) decompose EXACTLY into 4x4 rectangles
odd(Q_A) x odd(Q_B) with Q_A, Q_B 3-dim halves of class spaces (15 classes: 14 port classes + the normal class)."""
import json, sys
from collections import Counter, defaultdict
sys.path.insert(0,'../expBV'); from wordmodel import basis
P=json.load(open('instA.json'))['ports']; cen=json.load(open('cen.json'))['cen']; n=112
nor=[1,7,25,31,121,127,97,103]
dot=lambda a,b:(a&b).bit_count()&1
V=P+nor                       # 0..111 ports, 112..119 centre normals
cl=defaultdict(list)
for s in range(n): cl[tuple(cen[s])].append(s)
C=list(cl.values())+[list(range(112,120))]
cls={x:i for i,c in enumerate(C) for x in c}
ok=lambda c,m: print(('PASS ' if c else 'FAIL ')+m) or c
def span(B):
    s=[0]
    for v in B: s=s+[x^v for x in s]
    return s
ok(all(len(basis([V[x] for x in c]))==4 and sorted(V[x] for x in c)==sorted(y for y in span(basis([V[x] for x in c])) if y.bit_count()&1) for c in C),
   '15 classes (incl. normal class) = the 8 odd vectors of a 4-dim space')
rects=[]   # ordered (A,B): A tuple of 4 port idx (source side), B tuple of 4 idx in V (target side; >=112 = centre)
for i in range(14):
    for j in range(15):
        M={(a,b) for a in C[i] for b in C[j] if dot(V[a],V[b])==0}
        if not M: continue
        # bicliques: group rows by their neighbourhood
        rows=defaultdict(list)
        for a in C[i]: rows[tuple(sorted(b for b in C[j] if (a,b) in M))].append(a)
        for nb,As in rows.items():
            if not nb: continue
            rects.append((tuple(sorted(As)),tuple(nb)))
# centre rectangles: a whole port class (8) x 4 normals -> split along a split of that class already used elsewhere
used={A for A,B in rects if len(A)==4}
fix=[]
for A,B in rects:
    if len(A)==8:
        sp=[Q for Q in used if set(Q)<=set(A)]
        Q=sp[0]; Q2=tuple(sorted(set(A)-set(Q))); assert Q2 in used
        fix+= [(Q,B),(Q2,B)]
    else: fix.append((A,B))
rects=fix
pairs=Counter((a,b) for A,B in rects for a in A for b in B)
need={(s,t) for s in range(n) for t in range(n) if dot(P[s],P[t])==0}|{(s,112+q) for s in range(n) for q in range(8) if dot(P[s],nor[q])==0}
ok(set(pairs)==need and max(pairs.values())==1 and len(need)==5824+448, f'{len(rects)} ordered 4x4 rectangles partition the {len(need)} ordered pairs (5824 port + 448 centre)')
ok(all(len(A)==4==len(B) for A,B in rects),'all rectangles are 4x4')
isquad=lambda Q: len(basis([V[x] for x in Q]))==3
ok(all(isquad(A) and isquad(B) for A,B in rects),'every rectangle side = odd vectors of a 3-dim space (a quad node)')
ok(all(cls[A[0]]==cls[a] for A,B in rects for a in A) and all(cls[B[0]]==cls[b] for A,B in rects for b in B),'quads are class-pure (halves of class spaces)')
ga=Counter(A for A,B in rects); sc=Counter(B for A,B in rects)
print('gather quads used',len(ga),'out-degree dist',sorted(Counter(ga.values()).items()))
print('scatter quads used',len(sc),'in-degree dist',sorted(Counter(sc.values()).items()))
print('per class #gather quads',Counter(cls[A[0]] for A in ga),'\nper class #scatter quads',Counter(cls[B[0]] for B in sc))
json.dump(dict(V=V,C=C,rects=rects),open('quads.json','w'))
