"""expXP spectral Boot bound for PG(2,q) (numpy; about 10 s).
Theorem (field-free): cost(sigma) = 1/2 * total ascent of the cut profile >= 1/2 * cut(middle prefix) >= lambda_2 * n / 8 for the
incidence graph (n = 2(q^2+q+1) vertices, k = (q+1)/2), and lambda_2(PG(2,q)) = q + 1 - sqrt(q) (adjacency spectrum
+-(q+1), +-sqrt(q)).  So  Boot(PG(2,q)) >= (q+1-sqrt q)(q^2+q+1)/4,  Boot/E >= (1 - sqrt(q)/(q+1))/4.
Checks: (1) lambda_2 numerically on PG(2,p) for primes p <= 23 (homogeneous coordinates over F_p);
(2) the closed form against expXP's boot_pg.log; (3) the window statements over all prime powers q <= 10^4:
  LB > E/6 iff q >= 7 (PG(2,7): ceil(76.3) = 77 > 76), LB/E > thr(1/9) = 1/(1+5^(8/9)) = 0.19300 iff q >= 19
  (PG(2,17) misses by 1.5 sessions: 1065.05 < 1066.54).  Upper bounds quoted in the README are expXE/expXP greedy values.
Prints PASS_EXPXP_SPECTRAL."""
import itertools, math, numpy as np
ok = True
def chk(name, cond, info=''):
    global ok; ok &= bool(cond); print(('PASS ' if cond else 'FAIL ') + name, info, flush=True)

def pg_points(p):
    pts = []
    for v in itertools.product(range(p), repeat=3):
        if any(v):
            f = next(x for x in v if x)            # normalise: first nonzero coordinate = 1
            inv = pow(f, p - 2, p); pts.append(tuple(x * inv % p for x in v))
    return sorted(set(pts))

for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
    P = pg_points(p); N = len(P); idx = {x: i for i, x in enumerate(P)}
    A = np.zeros((2 * N, 2 * N))
    for i, x in enumerate(P):                     # lines = same normalised triples (duality)
        for j, y in enumerate(P):
            if sum(a * b for a, b in zip(x, y)) % p == 0: A[i, N + j] = A[N + j, i] = 1
    L = np.diag(A.sum(1)) - A; w = np.linalg.eigvalsh(L)
    lam2 = w[1]; n = 2 * N; E = (p + 1) * N; LB = lam2 * n / 8
    chk(f'PG(2,{p}): regular degree q+1, lambda_2 = q+1-sqrt(q)', set(A.sum(1)) == {p + 1} and abs(lam2 - (p + 1 - math.sqrt(p))) < 1e-9,
        f'n={n} E={E} lambda_2={lam2:.6f} LB={LB:.2f} LB/E={LB/E:.4f}')

LBf = lambda q: (q + 1 - math.sqrt(q)) * (q * q + q + 1) / 4
Ef = lambda q: (q + 1) * (q * q + q + 1)
log = {7: 76.3, 9: 159.25, 11: 288.72, 13: 475.55, 17: 1065.05, 19: 1489.81, 23: 2654.98}   # expXP boot_pg.log spec_lb
chk('closed form = boot_pg.log spectral LBs', all(abs(LBf(q) - v) < 0.01 for q, v in log.items()))
thr = 1 / (1 + 5 ** (8 / 9))
pp = [q for q in range(2, 10001) if (lambda f: len(f) == 1)({d for d in range(2, q + 1) if q % d == 0 and all(d % e for e in range(2, int(d ** .5) + 1))})]
above6 = [q for q in pp if math.ceil(LBf(q) - 1e-9) > Ef(q) / 6]
abovethr = [q for q in pp if LBf(q) > thr * Ef(q)]
chk('Boot > E/6 certified exactly for prime powers q >= 7', above6 == [q for q in pp if q >= 7], f'first: {above6[:4]}  PG(2,5): LB {LBf(5):.2f} -> {math.ceil(LBf(5))} vs E/6 = {Ef(5)/6:.0f}')
chk('Boot/E > 0.19300 (a = 1/9) certified exactly for prime powers q >= 19', abovethr == [q for q in pp if q >= 19],
    f'thr={thr:.5f}  q=17: LB {LBf(17):.2f} < {thr*Ef(17):.2f};  q=19: LB/E = {LBf(19)/Ef(19):.4f}')
chk('PG(2,7) window: ceil(LB) = 77 > E/6 = 76', math.ceil(LBf(7)) == 77 and Ef(7) / 6 == 76)
print('PASS_EXPXP_SPECTRAL' if ok else 'FAIL_EXPXP_SPECTRAL')
