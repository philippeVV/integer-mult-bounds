import sys; sys.path.insert(0,'.')
from sa_class import *
from collections import Counter
KA=[tuple(Ks[a]) for a in range(n)]
cnt=Counter(KA); print('distinct support words',len(cnt),'sizes',set(cnt.values()),'weights',set(len(k) for k in KA))
ws=sorted(cnt); print('pairwise intersections even:',all(len(set(x)&set(y))%2==0 for x in ws for y in ws))
# within a class: pairs even?
same=[(a,b) for a in range(n) for b in range(a+1,n) if KA[a]==KA[b]]
print('same-word pairs',len(same),'all even (q.q=1):',all(dot(P[a],P[b])==1 for a,b in same))
# odd pairs: overlap size distribution
print('overlap for odd pairs',Counter(len(set(KA[a])&set(KA[b])) for a in range(n) for b in range(n) if a<b and dot(P[a],P[b])==0))
print('overlap for even pairs',Counter(len(set(KA[a])&set(KA[b])) for a in range(n) for b in range(n) if a<b and dot(P[a],P[b])==1))
