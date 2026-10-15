"""F8 structure on F2^9 compatible with the dot product (self-dual basis), port F8-lines, free mask."""
import itertools, random
from common import *
MOD = 0b1011
def mul(a, b):
    r = 0
    for i in range(3):
        if b >> i & 1: r ^= a << i
    for i in (4, 3):
        if r >> i & 1: r ^= MOD << (i - 3)
    return r
def tr(x): return (x ^ mul(x, x) ^ mul(mul(x, x), mul(x, x))) & 7
SD = next(b for b in itertools.combinations(range(1, 8), 3)
          if all(tr(mul(b[i], b[j])) == (i == j) for i in range(3) for j in range(3)) and rank_f2(list(b)) == 3)
def h(u, v):
    s = 0
    for a, b in zip(u, v): s ^= mul(a, b)
    return s
def build(seed, trial):
    """Same sampling as b_graphs.py: trial 0 = identity structure."""
    random.seed(seed)
    for t in range(trial + 1):
        perm = list(range(9)); random.shuffle(perm)
        vs = [] if t == 0 else [random.choice([v for v in range(1, 512) if bin(v).count('1') % 2 == 0]) for _ in range(6)]
    if trial == 0: perm = list(range(9))
    P = load_ports()
    def mix(x):
        for v in vs:
            if dot(x, v): x ^= v
        return x
    def phi(x):
        v = [0, 0, 0]
        for k in range(9):
            if x >> k & 1: v[perm[k] // 3] ^= SD[perm[k] % 3]
        return v
    F = [phi(mix(q)) for q in P]
    n = len(P)
    free = [[a != b and h(F[a], F[b]) == 0 for b in range(n)] for a in range(n)]
    return F, free
