"""S(q) = dim W of the Hermitian Tanner code on the FULL plane PG(2,q), q = p^2 odd (expQK's family, odd char).

Point code: values f(phi_l) of degree-(k-1) forms, k = (q+1)/2, in coordinates of the lines through P
([q+1,k] doubly-extended RS, MDS).  Line code: same with Frobenius-twisted coordinates (x -> x^sqrt q).
At q = 9 this is the Braverman-He code (bh_build.py: degree-4 forms, cube twist at lines).
Singer/trace model: points g^i, lines phi_l(x) = Tr(g^-l x), edge (i,d), l = i+d, d in D = {d: Tr(g^-d) = 0}.
Block decomposition over characters of Z_N (with wrap cocycle rho, see README Sec.4):
    S = sum_chi [ (q+1) - rank [ G | -diag(rho^d w^{chi d}) H ] ].
usage: odd_singer.py p [dense]     (p = 3, 5, 7 -> q = 9, 25, 49)"""
import sys, time, math, collections
import numpy as np

p = int(sys.argv[1]); DENSE = len(sys.argv) > 2
q = p * p; e6 = 6; NF = p ** 6; ORD = NF - 1; N = q * q + q + 1; k = (q + 1) // 2; sig = p
assert ORD == N * (q - 1)
DIG = np.array([[(a // p ** t) % p for t in range(e6)] for a in range(NF)], dtype=np.int64)
PW = p ** np.arange(e6, dtype=np.int64)


def add(a, b): return ((DIG[a] + DIG[b]) % p) @ PW
def neg(a): return ((p - DIG[a]) % p) @ PW


def _pmulmod(a, b, tail):
    r = [0] * (2 * e6 - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): r[i + j] = (r[i + j] + x * y) % p
    for t in range(2 * e6 - 2, e6 - 1, -1):          # x^6 = -sum tail[s] x^s
        c = r[t]
        if c:
            r[t] = 0
            for s in range(e6): r[t - e6 + s] = (r[t - e6 + s] - c * tail[s]) % p
    return r[:e6]


def _xpow(n, tail):
    res = [1] + [0] * (e6 - 1); base = [0, 1] + [0] * (e6 - 2)
    while n:
        if n & 1: res = _pmulmod(res, base, tail)
        base = _pmulmod(base, base, tail); n >>= 1
    return res


def _primes(n):
    out, d = [], 2
    while d * d <= n:
        while n % d == 0: out.append(d); n //= d
        d += 1
    return sorted(set(out + ([n] if n > 1 else [])))


def find_prim():
    import itertools
    one = [1] + [0] * (e6 - 1)
    for tail in itertools.product(range(p), repeat=e6):
        if tail[0] == 0: continue
        if _xpow(ORD, list(tail)) != one or any(_xpow(ORD // r, list(tail)) == one for r in _primes(ORD)):
            continue
        # monic x^6 + sum tail[t] x^t ; build exp table by multiplying by x
        EXP = np.zeros(ORD, dtype=np.int64); v = np.zeros(e6, dtype=np.int64); v[0] = 1
        seen_one = False; ok = True
        for i in range(ORD):
            EXP[i] = int(v @ PW)
            if i > 0 and EXP[i] == 1: ok = False; break
            top = v[-1]; v = np.roll(v, 1); v[0] = 0
            v = (v - top * np.array(tail)) % p
        if ok and int(v @ PW) == 1: return EXP
    raise RuntimeError


EXP = find_prim(); LOG = np.full(NF, -1, dtype=np.int64); LOG[EXP] = np.arange(ORD)
assert (LOG[1:] >= 0).all()
EXP2 = np.concatenate([EXP, EXP])


def mul(a, b):
    a = np.asarray(a); b = np.asarray(b)
    r = EXP2[(LOG[a] + LOG[b]) % ORD]
    return np.where((a == 0) | (b == 0), 0, r)


def gp(n): return int(EXP[n % ORD])
def inv(a): return int(EXP[(-LOG[a]) % ORD])


def tr(a):
    s = 0
    for t in (1, q, q * q):
        s = int(add(s, int(EXP[(LOG[a] * t) % ORD])) if a else 0)
    return s


D = [d for d in range(N) if tr(gp(-d)) == 0]
assert len(D) == q + 1, len(D)
Fq = set(int(EXP[(N * t) % ORD]) for t in range(q - 1)) | {0}
d1, d2 = D[0], D[1]


def coords(b1, b2, y):
    """lam, mu in F_q with y = lam*b1 + mu*b2 (elements of F_{q^3}, F_q-linear)."""
    for lam in Fq:
        r = add(y, neg(int(mul(lam, b1))))
        # r = mu*b2 -> mu = r / b2 must lie in F_q
        mu = 0 if r == 0 else int(mul(r, inv(b2)))
        if mu in Fq: return lam, mu
    raise AssertionError


def fpow(a, n): return 1 if n == 0 else (0 if a == 0 else gp(LOG[a] * n))


def mono(lam, mu): return [int(mul(fpow(lam, k - 1 - r), fpow(mu, r))) for r in range(k)]


# point side: phi_d = lam phi_{d1} + mu phi_{d2}  <=>  g^-d = lam g^-d1 + mu g^-d2  (Tr nondegenerate)
G = np.array([mono(*coords(gp(-d1), gp(-d2), gp(-d))) for d in D], dtype=np.int64)
# line side: point g^-d on kernel(Tr) in basis (g^-d1, g^-d2), twisted coordinates
def twist(lm): return (fpow(lm[0], sig), fpow(lm[1], sig))
H = np.array([mono(*twist(coords(gp(-d1), gp(-d2), gp(-d)))) for d in D], dtype=np.int64)
theta_star_exp = N * (sig - 1) * (k - 1)            # theta* = g^{N (sigma-1)(k-1)}
rho_exp = (sig - 1) * (k - 1)                         # rho^N = theta*
omega_exp = q - 1                                     # omega = g^{q-1}, order N


def rank(M):
    M = M.copy(); r = 0; nr, nc = M.shape
    for c in range(nc):
        if r >= nr: break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0: continue
        i = r + nz[0]
        if i != r: M[[r, i]] = M[[i, r]]
        M[r] = mul(inv(int(M[r, c])), M[r])
        col = M[:, c].copy(); col[r] = 0; rows = np.nonzero(col)[0]
        if len(rows):
            sub = mul(col[rows][:, None], M[r][None, :])
            M[rows] = ((DIG[M[rows]] - DIG[sub]) % p) @ PW
        r += 1
    return r


def block_S():
    spec = []
    negH = ((p - DIG[H]) % p) @ PW
    Dn = np.array(D, dtype=np.int64)
    for chi in range(N):
        lam = EXP[((rho_exp + omega_exp * chi) * Dn) % ORD]
        M = np.hstack([G, mul(lam[:, None], negH)])
        spec.append((q + 1) - rank(M))
    return spec


def dense_S():
    """direct: unknowns a_i (k each), b_l (k each); equations G_d a_i = c(i,d) H_d b_l per edge."""
    E = N * (q + 1); M = np.zeros((E, 2 * N * k), dtype=np.int64); row = 0
    th = gp(theta_star_exp)
    for i in range(N):
        for di, d in enumerate(D):
            l = (i + d) % N; c = th if i + d >= N else 1
            M[row, i * k:(i + 1) * k] = G[di]
            hv = mul(c, H[di]); M[row, N * k + l * k:N * k + (l + 1) * k] = ((p - DIG[hv]) % p) @ PW
            row += 1
    return 2 * N * k - rank(M)


t0 = time.time()
print(f"q={q} p={p} N={N} k={k} |D|={len(D)} E={N*(q+1)}", flush=True)
sp = block_S(); S = sum(sp); E = N * (q + 1)
r = S / E; mk = (E - S) / S
a = 1 - math.log(mk) / math.log(5)
print(f"q={q}: S=dimW={S} spectrum={dict(sorted(collections.Counter(sp).items()))} E={E} S/E={r:.4f} "
      f"(E-S)/S={mk:.4f} a(b=5)={a:+.4f}  [{'>= 1/6' if r >= 1/6 else '< 1/6'}]  ({time.time()-t0:.1f}s)", flush=True)
if DENSE:
    Sd = dense_S(); print(f"dense check: S={Sd}  ({'agrees' if Sd == S else 'MISMATCH'})  ({time.time()-t0:.1f}s)", flush=True)

# limits-1 replay: closed form S(q) = [(q+1)(q^2+q+1) - sqrt(q)(q^2+1)]/4 (fitted on q = 9, 25; predicted q = 49)
assert 4 * S == (q + 1) * (q * q + q + 1) - p * (q * q + 1), S
assert not DENSE or Sd == S
print(f'PASS_EXPXJ_SQ{q}')
