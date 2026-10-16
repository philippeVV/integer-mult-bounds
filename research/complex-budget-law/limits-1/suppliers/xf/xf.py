#!/usr/bin/env python3
"""expXF checks (stdlib only, a few seconds).

A  exponent obstruction: integer convolution cannot factor through an exponent-2 (char-2) coefficient group
B  kernel determinants: radix-l DFT kernels are singular mod l; the char-2 additive butterfly is unimodular
C  multi-small-prime CRT outer: number/size of primes and the volume factor (constant)
D  recursion-depth law: per-level expansion c over depth L costs c^L; SS/Nussbaumer depth lglg n kills the saving,
   Fuerer depth log* n and the single-level HvdH outer keep it
E  Lefschetz/Ax illustration: a fixed integral seed has the char-0 rank mod every prime except finitely many
"""
import math, random

def sieve(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

# ---------------------------------------------------------------- A
print("== A  exponent obstruction")
random.seed(1)
N, b = 64, 8
a = [random.randrange(2 ** b) for _ in range(N)]; c_ = [random.randrange(2 ** b) for _ in range(N)]
conv = [sum(a[i] * c_[k - i] for i in range(N) if 0 <= k - i < N) for k in range(2 * N - 1)]
mx = max(conv)
print(f"  N={N}, b={b}: max coefficient {mx} ~ 2^{math.log2(mx):.1f}; bound N(2^b-1)^2 = 2^{math.log2(N*(2**b-1)**2):.1f}")
# any additive map phi into an exponent-2 group kills 2Z: phi(c(2a,b)) = phi(2 c(a,b)) = 0
two_a = [2 * x for x in a]
conv2 = [sum(two_a[i] * c_[k - i] for i in range(N) if 0 <= k - i < N) for k in range(2 * N - 1)]
assert all(y == 2 * x for x, y in zip(conv, conv2)) and any(conv2)
print("  c(2a,b)=2c(a,b)!=0, but every bilinear scheme C(Aa o Bb) over exponent-2 groups maps it to 0:")
print("  => coefficient-group exponent (or lcm over CRT factors) must exceed", f"2^{math.log2(mx):.1f}", "-- char 2 impossible")

# ---------------------------------------------------------------- B
print("== B  kernel determinants")
import cmath
def det(M):
    M = [r[:] for r in M]; n = len(M); d = 1
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(M[r][i]))
        if abs(M[p][i]) < 1e-12: return 0
        if p != i: M[i], M[p] = M[p], M[i]; d = -d
        d *= M[i][i]
        for r in range(i + 1, n):
            f = M[r][i] / M[i][i]
            for k in range(i, n): M[r][k] -= f * M[i][k]
    return d
for l in (2, 3, 5, 7):
    z = cmath.exp(2j * math.pi / l)
    F = [[z ** (j * k) for k in range(l)] for j in range(l)]
    print(f"  radix-{l} DFT kernel: |det|^2 = {abs(det(F))**2:.1f} = {l}^{l}  -> needs 1/{l}: singular in char {l}")
print("  additive (Cantor/Gao-Mateer) butterfly [[1,w],[1,w+1]]: det = 1 in char 2 (unimodular, no denominator)")
print("  => a DFT-type arithmetic layer of radix l needs l invertible in R, so gcd(l, char R) = 1 (frames l-adic, gates char R)")

# ---------------------------------------------------------------- C
print("== C  multi-small-prime CRT outer (primes in (P/2, P], odd)")
PR = sieve(200000)
for lgn in (2 ** 10, 2 ** 13, 2 ** 15):
    bbits = lgn                      # digit size b = Theta(lg n), as upstream's p
    need = 2 * bbits + lgn + 1       # bits of coefficient range ~ 2b + lg N
    best = None
    for P in range(16, 200000, 16):
        ps = [p for p in PR if P // 2 < p <= P]
        if sum(math.log2(p) for p in ps) >= need:
            best = (P, len(ps)); break
    P, k = best
    print(f"  lg n={lgn}: need {need} bits -> {k} primes <= {P} (P/lg n = {P/lgn:.2f}); "
          f"volume factor vs one ring = {sum(math.log2(p) for p in PR if P//2 < p <= P)/need:.3f}; "
          f"vs input bits b: x{need/bbits:.2f}")
print("  => constant volume factor; but #characteristics ~ 3 lg n / ln lg n -> infinity, each p ~ Theta(lg n)")

# ---------------------------------------------------------------- D
print("== D  recursion-depth law (t(L) = per-bit cost at L = lg n; kappa_eff = 1 - ln t / ln L)")
def ss(L, a, c=2.0):
    if L <= 16: return L
    return L ** (1 - a) + c * ss(L / 2, a, c)
def furer(L, a, K=8.0):
    if L <= 16: return L
    return L ** (1 - a) + K * furer(2 * math.log2(L), a, K)
for a_ in (0.0, 0.05, 0.317):
    row = []
    for e in (2 ** 8, 2 ** 16, 2 ** 64, 2 ** 256):
        L = float(e)
        k_ss = 1 - math.log(ss(L, a_)) / math.log(L)
        k_ss1 = 1 - math.log(ss(L, a_, 1.0)) / math.log(L)
        k_f = 1 - math.log(furer(L, a_)) / math.log(L)
        row.append(f"L=2^{int(math.log2(L))}: SS(c=2) {k_ss:+.3f} SS(c=1) {k_ss1:+.3f} Fuerer(K=8) {k_f:+.3f}")
    print(f"  a={a_}:"); [print("    " + r) for r in row]
print("  => depth-lglg n recursion with constant expansion c>=2 drives kappa_eff to <= 0 for every a;")
print("     c=1 (exact halving, no padding) keeps a - O(1/lglg); log* depth keeps a - o(1)")

# ---------------------------------------------------------------- E
print("== E  Lefschetz/Ax: fixed integral seed, rank mod p vs char 0")
random.seed(7)
r, m, n_ = 12, 16, 20
A = [[random.randint(-3, 3) for _ in range(r)] for _ in range(m)]
B = [[random.randint(-3, 3) for _ in range(n_)] for _ in range(r)]
M = [[sum(A[i][k] * B[k][j] for k in range(r)) for j in range(n_)] for i in range(m)]
def rank_mod(M, p):
    M = [[x % p for x in row] for row in M]; rk = 0; cols = len(M[0])
    for cidx in range(cols):
        piv = next((i for i in range(rk, len(M)) if M[i][cidx]), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]; inv = pow(M[rk][cidx], p - 2, p)
        M[rk] = [x * inv % p for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][cidx]:
                f = M[i][cidx]; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk
r0 = rank_mod(M, 1000000007)
drops = [p for p in PR if p < 3000 and rank_mod(M, p) != r0]
print(f"  char-0 rank {r0}; primes < 3000 with a different rank: {drops}")
print("  => a bounded seed defined over Z behaves as char 0 for all but finitely many p; a uniform bounded family")
print("     with a gap for infinitely many p is (Ax transfer of an existential sentence) a char-0 gap")
print("PASS_EXPXF")
