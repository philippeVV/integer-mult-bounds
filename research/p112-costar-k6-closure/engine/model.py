"""Independent rebuild of the E8 / P112 masked co-star model.

Line order follows PR #1320 MODEL.json (doubled coordinates):
  56 integral lines: for i<j lexicographic, +2 at i, then +2 / -2 at j;
  64 half-integral lines (1,s2..s8), tail signs from product((1,-1),repeat=7),
  even number of minus signs.
Inner product G(r,s) = v_r . v_s / 4  (so G(r,r)=2, off-diagonal in {0,+-1}).
F0 = [0,1,26,27,44,45,54,55] (e1+-e2, e3+-e4, e5+-e6, e7+-e8) is deleted.
"""
import itertools
from functools import lru_cache

def build_lines():
    L = []
    for i in range(8):
        for j in range(i + 1, 8):
            for s in (2, -2):
                v = [0] * 8
                v[i] = 2
                v[j] = s
                L.append(tuple(v))
    for tail in itertools.product((1, -1), repeat=7):
        if sum(1 for t in tail if t < 0) % 2 == 0:
            L.append((1,) + tail)
    return L

L = build_lines()
N = len(L)
assert N == 120
GM = [[sum(x * y for x, y in zip(L[a], L[b])) // 4 for b in range(N)] for a in range(N)]
for a in range(N):
    for b in range(N):
        assert sum(x * y for x, y in zip(L[a], L[b])) % 4 == 0
        assert GM[a][b] == (2 if a == b else GM[a][b]) and GM[a][b] in ((2,) if a == b else (0, 1, -1))

F0 = [0, 1, 26, 27, 44, 45, 54, 55]
LIVE = [i for i in range(N) if i not in F0]
NL = len(LIVE)
assert NL == 112
POS = {r: k for k, r in enumerate(LIVE)}   # source index -> live position 0..111
FULL = (1 << NL) - 1

def orth(a, b):
    return a != b and GM[a][b] == 0

# co-star masks over live positions
COSTAR = []
for a in range(N):
    m = 0
    for k, r in enumerate(LIVE):
        if GM[a][r] == 0 or r == a:
            m |= 1 << k
    COSTAR.append(m)

# orthogonality adjacency over live positions
ONB = [0] * NL
for x in range(NL):
    for y in range(NL):
        if x != y and GM[LIVE[x]][LIVE[y]] == 0:
            ONB[x] |= 1 << y

def bits(m):
    out = []
    while m:
        b = m & -m
        out.append(b.bit_length() - 1)
        m ^= b
    return out

def popc(m):
    return bin(m).count("1")

# ---------------- F2 helpers ----------------
def f2rank(rows):
    piv = {}
    rank = 0
    for r in rows:
        while r:
            h = r.bit_length() - 1
            if h in piv:
                r ^= piv[h]
            else:
                piv[h] = r
                rank += 1
                break
    return rank

def pairing_rank(rows_src, cols_src):
    """F2 rank of B[r,s]=(G(r,s)+1) mod 2 on rows_src x cols_src (source indices)."""
    rows = []
    for r in rows_src:
        m = 0
        for k, s in enumerate(cols_src):
            if (GM[r][s] + 1) % 2:
                m |= 1 << k
        rows.append(m)
    return f2rank(rows)

# ---------------- masked entries ----------------
def masked(a, b):
    """even entry (diagonal or orthogonal) between source indices"""
    return a == b or GM[a][b] == 0

def target(a, b):
    return 2 * GM[a][b]   # 4 on diagonal, 0 on orthogonal

if __name__ == "__main__":
    print("lines", N, "live", NL)
    print("index 56", L[56], "index 87", L[87])
    # masked entry count on P112
    cnt = sum(1 for a in LIVE for b in LIVE if masked(a, b))
    print("masked ordered entries on P112:", cnt)
    sizes = {}
    for a in range(N):
        sizes.setdefault(popc(COSTAR[a]), []).append(a)
    print("co-star sizes:", {k: len(v) for k, v in sizes.items()})
    print("size-56 normals:", sizes.get(56))
    # pairing rank of every co-star
    pr = set()
    for a in range(N):
        S = [LIVE[k] for k in bits(COSTAR[a])]
        pr.add(pairing_rank(S, S))
    print("co-star pairing ranks:", pr)
    print("full P112 pairing rank:", pairing_rank(LIVE, LIVE))
    print("orth degree in P112:", sorted(set(popc(ONB[x]) for x in range(NL))))
