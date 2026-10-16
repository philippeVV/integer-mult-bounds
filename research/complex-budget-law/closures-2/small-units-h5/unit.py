"""Exact cross-field unit toolkit (expNA).

A unit (complex orientation): v ports, pairing B over F2 (v x v, B_TT = 1), decoder N over Q given as a list of
totals (r_j, a_j) (integer/Fraction vectors), N = sum_j r_j a_j^T.  Exactness: N_TT != 0 and N_TS != 0 => B_TS = 0.
Window h = rank_F2 B.  Cost = sum_j rank_F2 B[supp r_j, supp a_j]  (#446 Thm C floor).  x = cost / v.
Bit orientation: pairing B over Q (Gram), decoder over F2; cost = sum_j rank_Q B[R_j, A_j].
"""
from fractions import Fraction
import itertools, math, random

# ---------------- F2 linear algebra on python-int bitsets ----------------
def rank_f2_rows(rows):
    """rows: list of ints (bitsets). returns rank."""
    basis = []  # list of (pivotbit, vec)
    r = 0
    piv = {}
    for x in rows:
        while x:
            hb = x.bit_length() - 1
            if hb in piv:
                x ^= piv[hb]
            else:
                piv[hb] = x
                r += 1
                break
    return r

def f2_matrix_rows(B, R, A):
    """B: list of row bitsets (bit s set iff B[T][s] = 1). Returns rows of block B[R, A] as bitsets over A positions."""
    Amask = 0
    for s in A:
        Amask |= 1 << s
    return [B[t] & Amask for t in R]

def rank_f2_block(B, R, A):
    return rank_f2_rows(f2_matrix_rows(B, R, A))

def mat_to_bits(M):
    """M: list of lists of 0/1 ints -> row bitsets"""
    out = []
    for row in M:
        b = 0
        for j, e in enumerate(row):
            if e & 1:
                b |= 1 << j
        out.append(b)
    return out

# ---------------- exact Q rank (fraction-free Bareiss on ints / Fractions) ----------------
def rank_q(M):
    """exact rank over Q of a matrix of ints or Fractions"""
    if not M or not M[0]:
        return 0
    # scale rows to integers
    rows = []
    for row in M:
        den = 1
        for e in row:
            if isinstance(e, Fraction):
                den = den * e.denominator // math.gcd(den, e.denominator)
        rows.append([int(e * den) for e in row])
    m, n = len(rows), len(rows[0])
    A = [r[:] for r in rows]
    rank = 0
    prev = 1
    col = 0
    rowi = 0
    for col in range(n):
        if rowi >= m:
            break
        p = None
        for i in range(rowi, m):
            if A[i][col] != 0:
                p = i
                break
        if p is None:
            continue
        A[rowi], A[p] = A[p], A[rowi]
        pv = A[rowi][col]
        for i in range(rowi + 1, m):
            if A[i][col] != 0 or True:
                ai = A[i][col]
                Ai = A[i]
                Ar = A[rowi]
                for k in range(col + 1, n):
                    Ai[k] = (pv * Ai[k] - ai * Ar[k]) // prev
                Ai[col] = 0
        prev = pv
        rowi += 1
        rank += 1
    return rank

def rank_mod_p(M, p=(1 << 61) - 1):
    """rank mod a large prime (<= rank_Q); entries ints or Fractions"""
    rows = []
    for row in M:
        rr = []
        for e in row:
            if isinstance(e, Fraction):
                rr.append(e.numerator % p * pow(e.denominator % p, p - 2, p) % p)
            else:
                rr.append(int(e) % p)
        rows.append(rr)
    m = len(rows)
    n = len(rows[0]) if m else 0
    r = 0
    for c in range(n):
        piv = None
        for i in range(r, m):
            if rows[i][c]:
                piv = i
                break
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = pow(rows[r][c], p - 2, p)
        for i in range(m):
            if i != r and rows[i][c]:
                f = rows[i][c] * inv % p
                Ri = rows[i]; Rr = rows[r]
                for k in range(c, n):
                    Ri[k] = (Ri[k] - f * Rr[k]) % p
        r += 1
        if r == m:
            break
    return r

# ---------------- unit evaluation ----------------
def supp(vec):
    return [i for i, e in enumerate(vec) if e != 0]

def decoder_matrix(totals, v):
    N = [[Fraction(0)] * v for _ in range(v)]
    for r, a in totals:
        Rs = supp(r); As = supp(a)
        for t in Rs:
            rt = r[t]
            row = N[t]
            for s in As:
                row[s] += rt * a[s]
    return N

def check_complex(Bbits, totals, v):
    """exact check of N∘B = diag with N over Q, B over F2. returns (ok, msg)"""
    N = decoder_matrix(totals, v)
    for t in range(v):
        if N[t][t] == 0:
            return False, 'zero diagonal at %d' % t
        if not (Bbits[t] >> t) & 1:
            return False, 'B diag zero at %d' % t
        for s in range(v):
            if s != t and N[t][s] != 0 and (Bbits[t] >> s) & 1:
                return False, 'lossy correction %d<-%d' % (t, s)
    return True, 'ok'

def cost_complex(Bbits, totals):
    return [rank_f2_block(Bbits, supp(r), supp(a)) for r, a in totals]

def window_f2(Bbits):
    return rank_f2_rows(list(Bbits))

def kappa(h, x):
    import format_map
    X, s = format_map.best(h, x)
    return X / (1 + X), s

# ---------------- root systems ----------------
def e8_roots():
    """240 roots of E8 scaled by 2 (integer coords; norm 8 -> true norm 2). returns list of tuples"""
    R = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (2, -2):
            for sj in (2, -2):
                v = [0] * 8; v[i] = si; v[j] = sj; R.append(tuple(v))
    for signs in itertools.product((1, -1), repeat=8):
        if signs.count(-1) % 2 == 0:
            R.append(tuple(signs))
    return R

def lines(roots):
    """one representative per +-pair"""
    seen = set(); out = []
    for r in roots:
        if tuple(-x for x in r) in seen:
            continue
        seen.add(r); out.append(r)
    return out

def ip(a, b):
    return sum(x * y for x, y in zip(a, b))

def gram(vecs, scale=1):
    return [[ip(a, b) // scale for b in vecs] for a in vecs]
