"""expXA lor.py: Lorentzian 'triple-root' family in I_{n,1} and the k-subset Frankl-Wilson family.

LZ_n : roots r_T = (1; 1_T), T a 3-subset of [n], form diag(-1, 1^n). Norm 2, pairwise ips in {-1,0,1}.
LZD_n: LZ_n plus the A_{n-1} roots e_i - e_j  (V = C(n+1,3)).  n = 8 gives 84 of the 120 E8 lines.
Decoder: N = W/2 = (1/2)(-y_0 y_0^T + sum_i y_i y_i^T), y_c = coordinate column c (n+1 terms).
FW_{n,k}: k-subsets (k even), W = (|a cap b|), labels (W+J) mod 2 (= odd vectors (1_a,1) of F2^{n+1}),
  decoder N = P(W) entrywise with P(t) = prod_{j<k/2} (t-2j)/(k-2j), expanded in binomials C(t,r) =
  sum_{|R|=r} [R in a][R in b]  -> one rank-1 term per r-set R with nonzero coefficient.
Exact validity (Fractions) and exact F2 costs.  usage: lor.py [nmax]
"""
import sys, itertools
from fractions import Fraction
from math import comb
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from common import label_bits, check_decoder, cost, rank_f2_rows


def lz(n, diffs=False, check=True):
    vecs = []
    for T in itertools.combinations(range(1, n + 1), 3):
        v = [0] * (n + 1); v[0] = 1
        for i in T: v[i] = 1
        vecs.append(v)
    if diffs:
        for i, j in itertools.combinations(range(1, n + 1), 2):
            v = [0] * (n + 1); v[i] = 1; v[j] = -1
            vecs.append(v)
    eps = [-1] + [1] * n
    V = len(vecs)
    W = [[sum(eps[c] * x[c] * y[c] for c in range(n + 1)) for y in vecs] for x in vecs]
    assert all(W[a][a] == 2 for a in range(V))
    assert all(W[a][b] in (-1, 0, 1) for a in range(V) for b in range(V) if a != b)
    Bb = label_bits(W)
    h = rank_f2_rows(Bb)
    terms = []
    for c in range(n + 1):
        p = {a: vecs[a][c] for a in range(V) if vecs[a][c]}
        terms.append((Fraction(eps[c], 2), p, p))
    if check:
        ok, info = check_decoder(W, Bb, terms)
        assert ok, info
    C, det = cost(Bb, terms)
    return V, h, C, det


def fw(n, k, check=True):
    sets = list(itertools.combinations(range(n), k))
    V = len(sets)
    idx = {s: i for i, s in enumerate(sets)}
    W = [[len(set(a) & set(b)) for b in sets] for a in sets]
    Bb = label_bits(W)
    h = rank_f2_rows(Bb)
    # P(t) = prod_{j<k/2} (t-2j)/(k-2j); binomial coefficients c_r via finite differences
    def P(t):
        num = Fraction(1)
        for j in range(k // 2):
            num *= Fraction(t - 2 * j, k - 2 * j)
        return num
    deg = k // 2
    c = []
    for r in range(deg + 1):  # c_r = Delta^r P(0)
        c.append(sum((-1) ** (r - i) * comb(r, i) * P(i) for i in range(r + 1)))
    terms = []
    for r in range(deg + 1):
        if c[r] == 0:
            continue
        for R in itertools.combinations(range(n), r):
            p = {i: 1 for i, s in enumerate(sets) if set(R) <= set(s)}
            terms.append((c[r], p, p))
    if check:
        ok, info = check_decoder(W, Bb, terms)
        assert ok, info
    C, det = cost(Bb, terms)
    return V, h, C, det, c


if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    print("LZ_n (triple roots in I_{n,1}), decoder W/2 (n+1 coordinate terms); exact check + F2 costs")
    for n in range(6, nmax + 1):
        for d in (False, True):
            V, h, C, det = lz(n, d, check=(n <= 14))
            assert not d or C == n * (n + 1), (n, C)      # LZD: n+1 coordinate terms of cost n each
            print(f" n={n:2d} {'LZD' if d else 'LZ '} V={V:5d} h={h:2d} terms={len(det)} C={C:4d} gamma={C/V:.4f} "
                  f"coord-terms={n+1} V/(h*terms)={V/(h*(n+1)):.3f} per-term={sorted(set(det))}")
    print("FW_{n,k} (k-subsets, Frankl-Wilson P(W) decoder)")
    for n, k in [(8, 4), (10, 4), (12, 4)]:
        V, h, C, det, cr = fw(n, k)
        print(f" n={n} k={k} V={V} h={h} terms={len(det)} C={C} gamma={C/V:.4f} c_r={[str(x) for x in cr]}")
    print('PASS_EXPXA_LOR')   # exact decoder checks are asserted inside lz()/fw()
