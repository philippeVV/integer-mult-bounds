"""Self-contained E8 line model for expQJ (independent of expNZ/e8.py's json).
Lines: 120 root lines (doubled coords, <x,y> = x.y/4). Labels u_r in F2^9 = (coords of r in a
simple-root basis mod 2, 1). Pairing B_rs = (1 + <r,s>) mod 2 (diag 1), checked = label form.
Term cost beta(R,S) = rank_F2 B[R,S]."""
import itertools
from fractions import Fraction as Fr

def _roots():
    out = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = [0] * 8; v[i] = 2 * si; v[j] = 2 * sj; out.append(tuple(v))
    for s in itertools.product((1, -1), repeat=8):
        if s.count(-1) % 2 == 0: out.append(tuple(s))
    return out

def canon(v):
    for a in v:
        if a: return tuple(v) if a > 0 else tuple(-x for x in v)

ROOTS = _roots(); assert len(ROOTS) == 240
LINES = sorted({canon(r) for r in ROOTS}); NL = 120; assert len(LINES) == NL
def ip(x, y):
    s = sum(a * b for a, b in zip(x, y)); assert s % 4 == 0; return s // 4
IP = [[ip(a, b) for b in LINES] for a in LINES]
BROW = [sum(1 << b for b in range(NL) if (1 + IP[a][b]) % 2) for a in range(NL)]
ORTH = [sum(1 << b for b in range(NL) if IP[a][b] == 0) for a in range(NL)]
FULL = (1 << NL) - 1

SIMPLE = [(1, -1, -1, -1, -1, -1, -1, 1), (2, 2, 0, 0, 0, 0, 0, 0), (-2, 2, 0, 0, 0, 0, 0, 0)] + \
         [tuple([0] * k + [-2, 2] + [0] * (6 - k)) for k in range(1, 6)]
CARTAN = [[ip(a, b) for b in SIMPLE] for a in SIMPLE]

def solve(M, b):
    """solve M x = b over Q (M square list of lists)"""
    n = len(M); A = [[Fr(x) for x in M[i]] + [Fr(b[i])] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0); A[c], A[p] = A[p], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [A[i][n] for i in range(n)]

# coordinates of each line in the simple-root basis (integers)
_ST = [list(col) for col in zip(*SIMPLE)]  # 8x8, columns = simple roots
COORD = []
for l in LINES:
    c = solve(_ST, l); assert all(x.denominator == 1 for x in c); COORD.append([int(x) for x in c])
LAB = [sum(((COORD[i][j] % 2) << j) for j in range(8)) | (1 << 8) for i in range(NL)]
assert len(set(LAB)) == NL
FORM = [[CARTAN[i][j] % 2 for j in range(8)] + [0] for i in range(8)] + [[0] * 8 + [1]]
def form(u, v):
    return sum(FORM[i][j] * ((u >> i) & 1) * ((v >> j) & 1) for i in range(9) for j in range(9)) % 2
assert all(((BROW[a] >> b) & 1) == form(LAB[a], LAB[b]) for a in range(NL) for b in range(NL))

def bits(m):
    while m:
        low = m & -m; yield low.bit_length() - 1; m ^= low
def rank2(rows):
    basis = []
    for r in rows:
        for b in basis: r = min(r, r ^ b)
        if r: basis.append(r)
    return len(basis)
def beta(R, S):
    return rank2([BROW[a] & S for a in bits(R)])
def span(R):
    return rank2([LAB[a] for a in bits(R)])
def supp(v, L=FULL):
    """lines a in L with <r_a, v> != 0 (v any rational vector in the doubled-coordinate frame)"""
    return sum(1 << a for a in bits(L) if sum(x * y for x, y in zip(LINES[a], v)) != 0)

# fundamental weights in doubled coordinates frame: omega_i with <omega_i, alpha_j> = delta (true ip = dot/4)
# solve SIMPLE . w = 4 e_i  (dot with doubled coords /4)
OMEGA = [solve(SIMPLE, [4 * int(i == j) for j in range(8)]) for i in range(8)]
def class_vec(mask):
    return [sum(OMEGA[i][t] for i in range(8) if (mask >> i) & 1) for t in range(8)]
CLASSES = [(m, supp(class_vec(m))) for m in range(1, 256)]
