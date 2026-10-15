"""expBK exact cone census on a raw word (hpe audit prefix, exact Fractions).
usage: census.py RAW.json.gz OUT.json
For every helper still ending at FULL: min over the 511 hyperplanes y^perp of the number of cone offenders
(row: supp(row h of P^-1) \\ {h}; column: g != h with P[g][h] != 0; own: h's end-of-B frame not inside y^perp)."""
import sys, json, pickle
from collections import Counter
from common import *

word = load(sys.argv[1])
v, R = word['report']['v'], word['report']['R']; n = 2*v + R
comp = compensation(word)
prefix = word['precondition'] + comp + word['A'] + word['cut'] + word['B']
P, Pinv = A.state_maps(prefix, n)
start, touches, endB = frames_trace(word)
hf = set(map(int, word.get('helper_final', {})))
col = {}
for g in range(n):
    for k in P[g]:
        col.setdefault(k, set()).add(g)
# positive control: current low-ended helpers pass the cone test against their own final frame
ctrl = Counter()
for h in sorted(hf):
    U = basis(word['helper_final'][str(h)])
    cone = (set(Pinv[h]) | col.get(h, set())) - {h}
    ctrl[(len(U), all(inside(endB[g], U) for g in cone) and inside(endB[h], U))] += 1
print('control low-ended (dim, pass):', sorted(ctrl.items()))
full = [h for h in range(2*v, n) if h not in hf]
Y = list(range(1, 512))
res = {}
dist = Counter()
for h in full:
    rowr = set(Pinv[h]) - {h}
    colr = col.get(h, set()) - {h}
    cone = rowr | colr
    best = None
    for y in Y:
        own = 0 if in_yperp(endB[h], y) else 1
        bad = sorted(g for g in cone if not in_yperp(endB[g], y))
        key = len(bad) + own
        if best is None or key < best[0]:
            best = (key, y, own, bad)
    k, y, own, bad = best
    res[h] = dict(k=k, y=y, own=own, offenders=[dict(g=g, row=g in rowr, col=g in colr, pinv=str(Pinv[h].get(g, 0)),
                                                    p=str(P[g].get(h, 0)), endB_dim=len(endB[g])) for g in bad],
                  cone=len(cone), endB_h_dim=len(endB[h]))
    dist[min(k, 30)] += 1
print('full-ended helpers', len(full), 'min offenders (incl own frame):', sorted(dist.items()))
json.dump({str(h): r for h, r in res.items()}, open(sys.argv[2], 'w'), indent=0)
pickle.dump(dict(P=P, Pinv=Pinv, comp=comp, endB=endB, start=start, touches=touches), open(sys.argv[2] + '.pkl', 'wb'))
for h, r in sorted(res.items()):
    if r['k'] == 1:
        print(h, 'y', r['y'], 'own', r['own'], 'endB_h', r['endB_h_dim'], 'cone', r['cone'], r['offenders'])
