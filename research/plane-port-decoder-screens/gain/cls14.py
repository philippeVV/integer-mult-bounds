import sys, json, itertools
from sa_class import *
from collections import defaultdict
r=json.load(open(sys.argv[1]))['r']
cls=defaultdict(list)
for a in range(n): cls[tuple(Ks[a])].append(a)
for key,C in cls.items():
    vs=[w for a in C for w in (P[a],r[a])]
    dimV=rank2(vs); rk=form_rank(vs)
    # q-span
    qs=[P[a] for a in C]
    prs=[[pair_rank([P[a],r[a]],[P[b],r[b]]) for b in C] for a in C]
    print(key,'dim span',dimV,'form rank',rk,'q-span dim',rank2(qs),'q form rank',form_rank(qs),'pair ranks',sorted(set(x for row in prs for x in row)))
