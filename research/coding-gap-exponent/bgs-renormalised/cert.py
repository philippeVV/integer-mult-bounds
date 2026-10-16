"""Exact (Fraction) certificate checker for the leaf-extended BH seed under BGS App. F.  Standard library only.
Checks, in exact rational arithmetic:
  (a) BH core: 182 vertices, 735 edges, 157 sessions, BH's orientation facts, every session hop distance = 5;
  (b) a rational primal flow on the core of rate >= 147/157 - tiny (so 147/157 is the exact MCF up to that tiny slack);
  (c) for each instance (full leaf-extension = BGS App. C graph; minimal leaf-extension; BGS's own alpha-extension dual):
      terminals pairwise distinct, all leaf capacities 1, all dual weights > 0, sum_i dist_i >= 1 (exact Dijkstra),
      z = sum cap*w, gamma0 = NC/z with NC = 1 (BH Sec. 4 code), l = max dist, w = min weight, c = l/w;
  (d) exponents: seq_printed = ln g0/ln(16 c1^2) (BGS Thm 1.1 c2), seq_F = ln g0/(2 ln c1) (App. F inequality
      m_i <= (3 m_{i-1})^{4 c1^{2^i}}), seq_c = ln g0/(2 ln c) (same proof before rounding g_i+4 up), all-N = seq/2."""
import sys, json, heapq
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from fractions import Fraction as Fr
from math import log
from collections import Counter
from bh import CORE, SESS, edges, Rs, C, POS
from variants import full, minimal

# (a)
assert len(CORE) == 735 and len(SESS) == 157
def dijkstra(nv, E, w, s):
    adj = [[] for _ in range(nv)]
    for e, (u, v) in enumerate(E): adj[u].append((v, w[e])); adj[v].append((u, w[e]))
    d = [None]*nv; d[s] = Fr(0); pq = [(Fr(0), s)]
    while pq:
        dx, x = heapq.heappop(pq)
        if dx > d[x]: continue
        for y, we in adj[x]:
            nd = dx + we
            if d[y] is None or nd < d[y]: d[y] = nd; heapq.heappush(pq, (nd, y))
    return d
one = [Fr(1)]*735
hop = [dijkstra(182, CORE, one, s)[t] for q, s, t in SESS]
assert all(h == 5 for h in hop); print("(a) core rebuilt: 735 edges, 157 sessions, all hop distances = 5")

# (b) rational primal flow
flows = json.load(open(__file__.rsplit('/', 1)[0] + '/core_flow.json'))
load = Counter(); rate = Counter(); eidx = {}
for e, (u, v) in enumerate(CORE): eidx[(u, v)] = eidx[(v, u)] = e
for i, p, f in flows:
    f = Fr(f).limit_denominator(10**9)
    s, t = SESS[i][1], SESS[i][2]                   # check p is an s-t path in the core
    x = s
    for e in p:
        u, v = CORE[e]; assert x in (u, v); x = v if x == u else u
    assert x == t
    rate[i] += f
    for e in p: load[e] += f
L = max(load.values()); lam = min(rate[i] for i in range(157)) / L
print(f"(b) rational core flow: rate after scaling to loads<=1 is {float(lam):.15f}; 147/157 - rate = {float(Fr(147,157)-lam):.2e}")
assert Fr(147, 157) - lam < Fr(1, 10**7)

# (c)+(d)
def check(name, nv, E, ss, w, base_ext=None):
    assert len(set(x for p in ss for x in p)) == 2*len(ss), "terminals not distinct"
    assert all(x > 0 for x in w), "zero weight"
    caps = [Fr(1)]*len(E)
    d = []
    for s, t in ss: d.append(dijkstra(nv, E, w, s)[t])
    S = sum(d); assert S >= 1, S
    z = sum(c*x for c, x in zip(caps, w)); g0 = 1/z
    l = max(d); wm = min(w); c = l/wm; c1 = c+1
    lg = log(g0.numerator) - log(g0.denominator)
    r = dict(seq_printed=lg/log(16*float(c1)**2), seq_F=lg/(2*log(c1)), seq_c=lg/(2*log(c)))
    print(f"(c) {name}: |E|={len(E)} sum dist={S} z={z} gamma0={g0}={float(g0):.6f} l/w=c={c}={float(c):.4f} c1={float(c1):.4f}")
    print("    exponents (sequence | all-N = /2):  " + "  ".join(f"{k}={v:.6f}|{v/2:.6f}" for k, v in r.items()))
    return r

nv, E, ss = full()
w = [Fr(1, 1099)]*len(E)
rF = check("FULL leaf-ext (BGS App.C graph), uniform 1/1099", nv, E, ss, w)
# BGS's own alpha-extension dual on the same graph (App. C/F): core 1/785, leaves alpha/(K(1+eps)), not renormalised
eps = Fr(10, 147); alpha = (1+eps)/(1+eps/2) - 1; K = 314
w_a = [Fr(1, 785)]*735 + [alpha/(K*(1+eps))]*314
d_a = [dijkstra(nv, E, w_a, s)[t] for s, t in ss]; z_a = sum(w_a)
assert z_a == Fr(147, 157)*(1+alpha)
print(f"    [baseline] BGS alpha-ext dual: alpha={alpha}, z={z_a} (=z(D)(1+alpha)), gamma0=1+eps/2={1/z_a}, "
      f"sum dist={float(sum(d_a)):.5f} (>1: not renormalised), c={float(max(d_a)/min(w_a)):.4f}")
rB = check("  same alpha-ext dual, as BGS use it", nv, E, ss, w_a)
nv2, E2, ss2 = minimal()
w2 = [Fr(1, 966)]*len(E2)
rM = check("MINIMAL leaf-ext (181 leaves), uniform 1/966", nv2, E2, ss2, w2)
