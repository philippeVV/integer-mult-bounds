#!/usr/bin/env python3
"""expNU / census.py -- exact (mod p = 2^61-1) separability census of the best realised helper network (expSL absorb10)
for idea 8: two sibling copies on identical labels, sources per copy, helper registers SHARED (one form over (x_A, x_B)).

Necessary condition used (Lemma S, README section 3):
  target t_A must form alpha_t x_A from what it reads on its chain: copy-A source forms (pure A), copy-B source forms
  (pure B, the only cross-copy side information) and shared helper forms (mixed).  Modulo the source forms, every
  helper read is a vector of the residual space Abar_t (+) Bbar_t, and t_A needs (abar_t, 0), t_B needs (0, abar_t).
  Hence: if abar_t != 0 then t needs >= 2 helper reads with independent residuals (single copy: >= 1).
  Per node X (no side information): usable_A(X) + usable_B(X) <= k_X, the number of registers arriving at X.
Reported: per-target helper reads h_t and single-copy residual rank r_t, per-node slack k_X - rank_X.
usage: python3 -I -B census.py [inputs/absorb10.json.gz inputs/absorb10_net.json.gz]
"""
import sys, json, os, gzip
from collections import Counter, defaultdict
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from geo import coef, n, LINE      # noqa: E402

p = (1 << 61) - 1
md = lambda q: (q.numerator * pow(q.denominator, -1, p)) % p
regf = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'inputs', 'absorb10.json.gz')
netf = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'inputs', 'absorb10_net.json.gz')
A = json.load(gzip.open(regf, 'rt')); D = json.load(gzip.open(netf, 'rt'))
E = [(tuple(u), tuple(v), f, tag) for u, v, f, tag in D['E']]


def rank(vs):
    Bs = []
    for v in vs:
        v = list(v)
        for piv, row in Bs:
            if v[piv]:
                c = v[piv]; v = [(x - c * y) % p for x, y in zip(v, row)]
        i = next((k for k, x in enumerate(v) if x), None)
        if i is None:
            continue
        inv = pow(v[i], -1, p); v = [x * inv % p for x in v]
        Bs.append((i, v))
    return len(Bs), Bs


def reduce(Bs, v):
    v = list(v)
    for piv, row in Bs:
        if v[piv]:
            c = v[piv]; v = [(x - c * y) % p for x, y in zip(v, row)]
    return v


# ---- assign every register step to one edge of E (sources prefer an edge whose form contains x_s)
pool = defaultdict(list)
for i, (u, v, f, tag) in enumerate(E):
    pool[(u, v)].append(i)
role_of = {}
order = sorted(range(len(A['pieces'])), key=lambda r: 0 if isinstance(A['role'][r], list) and A['role'][r][0] == 'S' else 1)
unmatched = 0
for r in order:
    rl = A['role'][r]
    kind = rl[0] if isinstance(rl, list) else rl
    if kind == 'T':
        continue
    nodes = [tuple(x) for x in A['pieces'][r]]
    for u, v in zip(nodes, nodes[1:]):
        cand = pool.get((u, v), [])
        if not cand:
            unmatched += 1; continue
        pick = cand[0]
        if kind == 'S':
            s = rl[1]
            pick = next((i for i in cand if E[i][2][s]), cand[0])
        cand.remove(pick)
        role_of[pick] = kind
left = sum(len(c) for c in pool.values())
print(f'edges {len(E)}; register steps matched {len(role_of)}; unmatched steps {unmatched}; edges without a register {left}')
print('edge roles', Counter(role_of.values()))

# ---- per-node slack (all in-edges, roles mixed)
inn = defaultdict(list)
for i, (u, v, f, tag) in enumerate(E):
    inn[v].append(i)
slack = Counter(); tot_k = tot_r = 0
for X, ids in inn.items():
    k = len(ids); r, _ = rank([E[i][2] for i in ids])
    tot_k += k; tot_r += r; slack[k - r] += 1
print(f'per-node: in-edges {tot_k}, sum of in-span ranks {tot_r}, dependent (slack) edges {tot_k - tot_r}; nodes by slack {dict(sorted(slack.items()))}')

# ---- per-target separability census
hist_h = Counter(); hist_r = Counter(); need = 0; fail_h = 0; fail_r = 0
for t, ch in D['tchains'].items():
    t = int(t)
    S, H = [], []
    for N in ch:
        N = tuple(N)
        if len(N) == 1:                         # own line node never on a target chain; keep for safety
            pass
        for i in inn.get(N, []):
            (S if role_of.get(i) == 'S' else H).append(E[i][2])
    _, Bs = rank(S)
    abar = reduce(Bs, [md(coef(t, s)) for s in range(n)])
    if not any(abar):
        hist_h['abar=0'] += 1; continue
    need += 1
    hres = [reduce(Bs, h) for h in H]
    r, _ = rank(hres)
    hist_h[len(H)] += 1; hist_r[r] += 1
    fail_h += len(H) < 2
    fail_r += r < 2
print(f'targets with abar_t != 0 (need helper content): {need}/112')
print('helper reads h_t on chain (histogram)', dict(sorted(hist_h.items(), key=str)))
print('single-copy residual rank r_t of those reads (histogram)', dict(sorted(hist_r.items())))
print(f'targets violating the joint necessary condition h_t >= 2: {fail_h}; with single residual rank < 2: {fail_r}')
assert len(E) == len(role_of) and unmatched == 0 and need == 112 and fail_h == 0 and fail_r == 0
assert tot_k == tot_r                     # zero dependence slack at every node
print('PASS_EXPNU_CENSUS')
