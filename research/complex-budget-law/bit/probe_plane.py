#!/usr/bin/env python3
"""Step 0 of the proposed bit experiment: plane-port upgrades of P111 ports sharing one second direction r.

usage: python3 -I -B probe_plane.py [--root REPO_ROOT]
For a root line r and a set S of live ports with q_a . r = 0, give each a in S the plane E_a = <q_a, r> (source E_a -> FULL,
target 0 -> E_a^perp); other ports stay lines. A decoder entry (i, j), i != j, is free only if E_i is orthogonal to E_j.
The line decoder B0 = (1 + G) mod 2 (G = <x,y>/2) has B0_ij = 1 on orthogonal pairs, so every orthogonal pair that loses
plane orthogonality is a newly forced zero. Keeping B0's other choices, the additive correction C must be supplied by new
copy terms: #terms >= rank_F2(C), each term costs >= 1. Gain of upgrading S: (8/9)|S| deficit units per invocation.
Pairwise-odd S avoids S x S corrections (r.r != 0 makes all S x S pairs non-orthogonal).
Standard library only.
"""
import gzip, hashlib, json, sys
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = Path(sys.argv[sys.argv.index('--root') + 1]) if '--root' in sys.argv else HERE.parents[2]
WORD = ROOT / 'research/e8-bit-p111-tuned/tune/bitE8mut-p111-tuned.json.gz'
assert hashlib.sha256(WORD.read_bytes()).hexdigest() == 'e3e301e62bc2b41dbf1ceac0acf3d35824ebeab1848c5a2142d4577b109da52f'
w = json.load(gzip.open(WORD, 'rt'))
live = [tuple(x) for x in w['lines']]
n = len(live)
dot = lambda a, b: sum(x * y for x, y in zip(a, b))
assert n == 111 and {dot(q, q) for q in live} == {4}

roots = set(live) | {tuple(-x for x in q) for q in live}
grown = True
while grown:
    grown = False
    for q in list(roots):
        for x in list(roots):
            y = tuple(int(xi - Q(2 * dot(x, q), dot(q, q)) * qi) for xi, qi in zip(x, q))
            if y not in roots:
                roots.add(y); grown = True
canon = lambda v: v if next(x for x in v if x) > 0 else tuple(-x for x in v)
lines = sorted({canon(r) for r in roots})
assert len(roots) == 240 and len(lines) == 120 and len(set(lines) - {canon(q) for q in live}) == 9


def rank2(vs):
    piv, r = {}, 0
    for v in vs:
        while v:
            h = v.bit_length() - 1
            if h in piv:
                v ^= piv[h]
            else:
                piv[h] = v; r += 1; break
    return r


B0 = [sum(((1 + dot(live[i], live[j]) // 2) % 2) << j for j in range(n)) for i in range(n)]
assert rank2(B0) == 9


def correction(S, r):
    """explicit forced-zero correction matrix (rows as bitmasks) for upgraded set S with second direction r."""
    E = {i: [live[i]] + ([r] if i in S else []) for i in range(n)}
    return [sum(1 << j for j in range(n) if j != i and dot(live[i], live[j]) == 0
                and any(dot(x, y) for x in E[i] for y in E[j])) for i in range(n)]


stats, sizes, best, checked = Counter(), Counter(), None, 0
for r in lines:
    perp = [a for a in range(n) if dot(live[a], r) == 0 and canon(live[a]) != r]
    T = {a: sum(1 << j for j in range(n) if j != a and dot(live[a], live[j]) == 0 and dot(live[j], r) != 0) for a in perp}
    for a in perp:
        sizes[bin(T[a]).count('1')] += 1
    adj = {a: {b for b in perp if b != a and dot(live[a], live[b]) != 0} for a in perp}
    cliques = []

    def bk(Rs, P, X):
        if not P and not X:
            cliques.append(frozenset(Rs)); return
        for v in list(P):
            bk(Rs | {v}, P & adj[v], X & adj[v]); P = P - {v}; X = X | {v}

    bk(set(), set(perp), set())
    for S in cliques:
        rk = 2 * rank2([T[a] for a in S])          # C = [[0, M], [M^T, 0]] with S and the T's disjoint
        if checked < 40:                            # direct check of the block identity on a sample
            assert rank2(correction(S, r)) == rk; checked += 1
        stats[(len(S), rk)] += 1
        g = Q(8, 9) * len(S) - rk
        best = g if best is None or g > best else best
out = dict(lines=120, live=n, B0_rank_F2=9,
           single_port_newly_forced_zeros_per_row=dict(sorted(sizes.items())),
           maximal_pairwise_odd_common_r_sets={f'|S|={s}, rank_F2(C)={k}': c for (s, k), c in sorted(stats.items())},
           max_gain_minus_rank=str(best), sampled_direct_checks=checked,
           verdict='additive repair of a common-r plane upgrade never pays: rank_F2(C) >= 1.6|S| > (8/9)|S|')
assert best < 0
print(json.dumps(out, indent=1))
print('PASS_BIT_PLANE_UPGRADE_STEP0')
