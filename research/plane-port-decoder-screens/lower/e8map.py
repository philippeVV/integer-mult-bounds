"""Exact identification of the 112 ports (+8 missing frame vectors) with E8 root lines:
find real roots rho_a (norm 2) with q_a.q_b = <rho_a,rho_b> + 1 (mod 2). Writes e8map.json."""
import itertools, json
from fractions import Fraction as Fr
from ports import P, dot
TARGET = [x for x in range(512) if bin(x).count('1') in (3, 7)]
assert len(TARGET) == 120 and set(P) <= set(TARGET)
MISSING = sorted(set(TARGET) - set(P))
# E8 roots (doubled coordinates, so inner product /4)
roots = []
for i, j in itertools.combinations(range(8), 2):
    for si in (2, -2):
        for sj in (2, -2):
            v = [0] * 8; v[i] = si; v[j] = sj; roots.append(tuple(v))
for signs in itertools.product((1, -1), repeat=8):
    if signs.count(-1) % 2 == 0: roots.append(signs)
ip = lambda a, b: sum(x * y for x, y in zip(a, b)) // 4
lines = []
for r in roots:
    if tuple(-x for x in r) not in lines: lines.append(r)
assert len(lines) == 120
# odd-ip graph on lines vs dot==0 graph on TARGET: backtracking isomorphism respecting linearity is
# replaced by plain graph isomorphism with refinement (both are the same 120-vertex graph).
A1 = {i: {j for j in range(120) if j != i and ip(lines[i], lines[j]) % 2} for i in range(120)}
A2 = {i: {j for j in range(120) if j != i and dot(TARGET[i], TARGET[j]) == 0} for i in range(120)}
assert all(len(A1[i]) == 56 for i in A1) and all(len(A2[i]) == 56 for i in A2)
def iso():
    order = []
    # BFS order on graph 1 for good pruning
    seen = [0]; order = [0]
    while len(order) < 120:
        best = max((j for j in range(120) if j not in seen), key=lambda j: len(A1[j] & set(order)))
        seen.append(best); order.append(best)
    m = {}
    used = set()
    def bt(k):
        if k == 120: return True
        i = order[k]
        for c in range(120):
            if c in used: continue
            if all((m[j] in A2[c]) == (j in A1[i]) for j in order[:k]):
                m[i] = c; used.add(c)
                if bt(k + 1): return True
                del m[i]; used.discard(c)
        return False
    assert bt(0)
    return m
m = iso()
# check: the frame (pairwise even ip, i.e. orthogonal lines) set mapping to MISSING may need an automorphism;
# instead report which lines map to MISSING and verify they form a frame (pairwise orthogonal roots).
inv = {TARGET[c]: lines[i] for i, c in m.items()}
fr = [inv[x] for x in MISSING]
assert all(ip(a, b) == 0 for a, b in itertools.combinations(fr, 2)), 'missing set is not a frame'
for x in TARGET:
    for y in TARGET:
        assert (dot(x, y) - ip(inv[x], inv[y]) - 1) % 2 == 0
json.dump({str(x): list(inv[x]) for x in TARGET}, open('e8map.json', 'w'))
print('OK e8 identification; frame of missing lines verified; all 14400 pair parities match')
