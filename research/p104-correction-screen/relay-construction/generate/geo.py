"""Load P104 geometry (same code path as check_p104_window_bill.py) and dump a compact JSON."""
import json, os, sys, collections
from fractions import Fraction as Q
from itertools import combinations, product
ROOT = sys.argv[1]
PKG = os.path.join(ROOT, 'research', 'p104-native-channel')
W = json.load(open(os.path.join(PKG, 'inputs/LOWRANK-WITNESS.json')))
fx = json.load(open(os.path.join(PKG, 'inputs/P104-FIXTURE.json')))['section4_explicit_P104']
L = set()
for i, j in combinations(range(8), 2):
    for sg in (-1, 1):
        v = [0]*8; v[i] = 2; v[j] = 2*sg; L.add(tuple(v))
for sg in product((-1, 1), repeat=7):
    if sg.count(-1) % 2 == 0: L.add((1, *sg))
removed = {tuple(r) for r in fx['F0'] + fx['C']}
ports = [r for r in sorted(L) if r not in removed]; n = len(ports)
def span(vs):
    s = {0}
    for b in vs:
        if b not in s: s |= {x ^ b for x in s}
    return frozenset(s)
cur = {r: span(f) for r, f in W['initial_frames'].items()}
path = {r: [(0, cur[r])] for r in cur}
for i, e in enumerate(W['program']):
    if e[0] == 'transport':
        assert cur[e[1]] == span(e[2]); new = span(e[3]); cur[e[1]] = new; path[e[1]].append((i + 1, new))
    elif e[0] == 'zero': cur['C'] = span(e[2]); path['C'] = [(i + 1, cur['C'])]
X = ['X%d' % s for s in range(n)]; Y = ['Y%d' % t for t in range(n)]
q = [W['initial_frames'][X[s]][0] for s in range(n)]
H = [cur[Y[t]] for t in range(n)]
def normal(Hs):
    ns = [v for v in range(1, 512) if all(bin(v & h).count('1') % 2 == 0 for h in Hs)]
    assert len(ns) == 1; return ns[0]
def basis(fs):
    b = []
    for v in sorted(fs):
        w = v
        for x in b: w = min(w, w ^ x)
        if w: b.append(w)
    return b
NE = len(W['program'])
def segs(r):
    p = path[r]; out = []
    for k, (c, f) in enumerate(p):
        end = p[k+1][0] if k + 1 < len(p) else NE + 1
        out.append(dict(start=c, end=end, basis=basis(f), dim=len(basis(f))))
    return out
out = dict(n=n, NE=NE, ports=ports, q=q, Hnormal=[normal(h) for h in H],
           G=[[sum(x*y for x, y in zip(a, b))//4 for b in ports] for a in ports],
           X=[segs(x) for x in X], Y=[segs(y) for y in Y],
           final_X=[basis(cur[x]) for x in X])
json.dump(out, open(sys.argv[2], 'w'))
print('ok', n, NE)
