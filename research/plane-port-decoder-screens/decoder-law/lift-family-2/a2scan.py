import collections
from a2lib import *
K = ylines_T()
print('|T0|', len(T0), '|T_i|', [len(Ti[x]) for x in S])
cnt = collections.Counter((len(k), tuple(sorted(x for x in CLASS[h] if x in S))) for h, k in K.items())
print('(kernel dim of Y[O cap T], O cap S) counts:', sorted(cnt.items()))
lines = collections.defaultdict(list)
for h, k in K.items():
    if len(k) == 1: lines[normline(k[0])].append(h)
print('distinct candidate lines', len(lines))
for g, hs in lines.items():
    print('  line', [str(x) for x in g], 'classes', [(h, tuple(x for x in CLASS[h] if x in S)) for h in hs])
