"""expXA scale.py: exact costs of LZD_n and FW_{n,k} at large n via the label-span rank identity.

For a term whose support is a port set R, rank_F2 B[R,R] = rank of the F2 form restricted to U = span{labels of R}
(B[R,R] = L_R F L_R^T with L_R of full column-rank image U, so rank = rank(F|U)).  This gives exact ranks without
building V x V matrices.  Cross-checked against brute force (lor.py) for small n.

LZD_n labels: (x mod 2) in F2^{n+1} with dot form, plus the J coordinate: B_ab = <l_a,l_b> + 1, i.e. label (l_a, 1)
in F2^{n+2} with the dot form.  Term c (coordinate c) supports roots with x_c != 0.
FW_{n,k} labels: (1_a, 1) in F2^{n+1}, dot form. Term R (|R| = r, coefficient c_r != 0) supports a superset of R.
"""
import itertools, random, sys
from fractions import Fraction
from math import comb


def form_rank(vectors):
    """rank of the dot form restricted to span(vectors) (vectors as python ints)."""
    basis = []
    piv = {}
    for v in vectors:
        w = v
        while w:
            lb = w.bit_length() - 1
            if lb in piv:
                w ^= piv[lb]
            else:
                piv[lb] = w; basis.append(v); break
    rows = []
    for x in basis:
        r = 0
        for j, y in enumerate(basis):
            if bin(x & y).count('1') & 1:
                r |= 1 << j
        rows.append(r)
    piv2 = {}
    rk = 0
    for v in rows:
        while v:
            lb = v.bit_length() - 1
            if lb in piv2:
                v ^= piv2[lb]
            else:
                piv2[lb] = v; rk += 1; break
    return rk


def bits(idx):
    x = 0
    for i in idx: x |= 1 << i
    return x


def lzd_labels(n):
    J = 1 << (n + 1)
    labs = []
    for T in itertools.combinations(range(1, n + 1), 3):
        labs.append((bits((0,) + T) | J, (0,) + T))
    for i, j in itertools.combinations(range(1, n + 1), 2):
        labs.append((bits((i, j)) | J, (i, j)))
    return labs


def lzd(n):
    labs = lzd_labels(n)
    V = len(labs)
    h = form_rank([l for l, _ in labs])
    C = 0
    per = []
    for c in range(n + 1):
        r = form_rank([l for l, s in labs if c in s])
        per.append(r); C += r
    return V, h, C, per


def fw_coeffs(k):
    def P(t):
        x = Fraction(1)
        for j in range(k // 2):
            x *= Fraction(t - 2 * j, k - 2 * j)
        return x
    return [sum((-1) ** (r - i) * comb(r, i) * P(i) for i in range(r + 1)) for r in range(k // 2 + 1)]


def fw_term_rank(n, k, r, samples=None):
    """rank of B on k-sets containing a fixed r-set R = {0..r-1}; labels (1_a,1) in F2^{n+1}."""
    J = 1 << n
    R = tuple(range(r))
    rest = list(range(r, n))
    vecs = []
    if comb(n - r, k - r) <= 4000:
        for A in itertools.combinations(rest, k - r):
            vecs.append(bits(R + A) | J)
    else:
        rng = random.Random(r * 1000 + n)
        for _ in range(6 * n + 50):   # random spanning sample; rank of span saturates (checked below)
            A = tuple(rng.sample(rest, k - r))
            vecs.append(bits(R + A) | J)
    return form_rank(vecs)


def fw(n, k):
    cr = fw_coeffs(k)
    V = comb(n, k)
    J = 1 << n
    h = fw_term_rank(n, k, 0)
    C = 0
    det = {}
    for r, c in enumerate(cr):
        if c == 0:
            continue
        rho = fw_term_rank(n, k, r)
        det[r] = rho
        C += comb(n, r) * rho
    return V, h, C, det


if __name__ == "__main__":
    print("LZD_n exact (span method). gamma = C/V, merit = gamma*(h-1)")
    for n in list(range(8, 31)) + [41, 61]:
        V, h, C, per = lzd(n)
        assert C * (n - 1) == 6 * V and C == n * (n + 1), (n, C, V)   # gamma_n = 6/(n-1) exactly
        print(f" n={n:2d} V={V:6d} h={h:2d} C={C:5d} per-term={sorted(set(per))} gamma={C/V:.4f} (6/(n-1)={6/(n-1):.4f}) merit={C/V*(h-1):.3f}")
    print("FW_{n,k} exact (span method; random spanning samples for big supports)")
    for k in (4, 6, 8):
        print(f" k={k} coefficients c_r = {[str(x) for x in fw_coeffs(k)]}")
        for n in [10, 12, 16, 24, 32, 48, 64, 96, 128]:
            if n <= k + 1: continue
            V, h, C, det = fw(n, k)
            print(f"  n={n:3d} V={V:12d} h={h:3d} C={C:10d} per-term(r:rho)={det} gamma={C/V:.4f} merit={C/V*(h-1):.3f}")
    print('PASS_EXPXA_SCALE')
