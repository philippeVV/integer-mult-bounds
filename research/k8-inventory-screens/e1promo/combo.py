"""Joint minimal completions of E1->E2 promotions: pairs, groups, all; greedy and exhaustive over small sets (expB).
usage: combo.py WORD cones analysis.json"""
import sys, json, itertools
from collections import Counter
sys.path.insert(0, __import__('os').path.dirname(__file__))
import e1lib as L
w, sha = L.load(sys.argv[1]); W = L.Word(w); W.sha = sha; W.prepare(); W.compute_cones(sys.argv[2]); C, gc = W.ceilings()
an = json.load(open(sys.argv[3]))
cands = {}  # helper -> list of E' masks
for q, r in an['helpers'].items():
    for c in r['candidates']:
        if c['legal']: cands.setdefault(int(q), []).append(L.mask(c['Ep']))
H = sorted(cands)
width = {q: L.dim(W.end[q]) - 1 for q in H}
def evalset(sel):
    ns = {q: m for q, m in sel.items()}
    raised, dprof, legal, bad, touched = W.complete(ns)
    dw = Counter()
    for q in sel: dw[width[q] - 1] += 1; dw[width[q]] -= 1
    drc, parts = L.price(dprof, dw)
    return drc, dprof, legal, raised, parts
print('promotable', len(H), 'candidates per helper', Counter(len(v) for v in cands.values()))
single = {q: min((evalset({q: m})[0], i) for i, m in enumerate(cands[q])) for q in H}
best1 = min(single.items(), key=lambda kv: kv[1][0]); print('best single', best1)
# all jointly (first-choice E')
allsel = {q: cands[q][single[q][1]] for q in H}
drc, dprof, legal, raised, parts = evalset(allsel)
print('ALL', len(allsel), 'legal', legal, 'dRC %.6f' % drc, dict(sorted(dprof.items())), parts, 'raised', len(raised))
# pairs: interactions
inter = []
for a, b in itertools.combinations(H, 2):
    d2 = evalset({a: cands[a][single[a][1]], b: cands[b][single[b][1]]})[0]
    gain = d2 - single[a][0] - single[b][0]
    if abs(gain) > 1e-12: inter.append((gain, a, b, d2))
inter.sort(); print('interacting pairs', len(inter), inter[:10])
# greedy
cur = {}; curv = 0.0
while True:
    best = None
    for q in H:
        if q in cur: continue
        for m in cands[q]:
            s = dict(cur); s[q] = m; d = evalset(s)[0]
            if best is None or d < best[0]: best = (d, q, m)
    if best is None or best[0] >= curv - 1e-12: break
    cur[best[1]] = best[2]; curv = best[0]; print('greedy add', best[1], 'dRC %.6f' % curv, flush=True)
print('greedy final', sorted(cur), 'dRC %.6f' % curv)
