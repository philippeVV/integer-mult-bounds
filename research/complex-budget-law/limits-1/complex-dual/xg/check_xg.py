"""expXG checker.  run:  uv run -q --no-project --python 3.12 --with numpy --with scipy python -I -B check_xg.py   (numpy+scipy; ~1 min; prints PASS_EXPXG)
C1 free values: flat Fourier codes give Theta_2 = 2/11 (PG3, girth 6) and Theta_3 = 1/8 (W3, girth 8).
C2 Theorem B (Barthe gauge + T1): every MDS half-rate code can be gauged so that one side is flat -> m1 = 1/4,
   hence dim W <= E/4 on every graph; tight on K_{4,4} (product code).
C3 Theorem C (bootstrap, any field): PG3 has dim W <= 8; expXE's r8 optimum is a genuine MDS code with dim W = 8 (tight).
C4 Lemma D (gauge artefact): for ANY codes, the gauge t^{1_{X^c}} (X a 2-factor) drives the |E|/2 smallest
   eigenvalues of K = 2I - Pi_P - Pi_L to 0 while dim W is unchanged; eigenvalue losses certify nothing.
C5 Theorem E (product): (U_P x U'_P) cap (U_L x U'_L) = W x W', Theta of the product = product of Thetas.
"""
import numpy as np, itertools, os, sys
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from xg_core import *
from boot import boot
ok = True
def check(cond, msg):
    global ok
    print(('ok   ' if cond else 'FAIL ') + msg); ok &= bool(cond)

# C1 (+ Theorem F): free (arcsine) moments; moment-only Markov = expRF's M_g, kernel-atom refinement = exactly 1/g
from math import comb
for J in range(1, 8):
    m = np.array([1.0] + [comb(2 * k, k) / 4 ** k / 2 for k in range(1, J + 1)])
    check(abs(markov_ref(m, J) - 1 / (2 * J + 2)) < 1e-9, f'C1 free law, J={J} (girth {2*J+2}): moment-only {markov(m, J):.6f}, refined {markov_ref(m, J):.6f} = 1/g')
V, nP, E, g = load('PG3'); Qf = fourier_codes(V, np.random.default_rng(0))
PPf, PLf = projectors(V, nP, E, Qf); mf = moments(PPf, PLf, 2)
check(abs(markov(mf, 2) - 2 / 11) < 1e-9 and abs(markov_ref(mf, 2) - 1 / 6) < 1e-9, 'C1 PG3 Fourier codes: Theta_2 = 2/11 (expRF), refined = 1/6 = threshold')
V8, nP8, E8, g8 = load('W3'); Qf8 = fourier_codes(V8, np.random.default_rng(0))
check(abs(theta(V8, nP8, E8, Qf8, 3) - 1 / 8) < 1e-9, 'C1 W3 Fourier Theta_3 = 1/8 = %.6f' % theta(V8, nP8, E8, Qf8, 3))

# C2: Barthe flat gauge on random non-flat MDS codes, several graphs
for name in ['PG3', 'BA4', 'W3']:
    Vx, nPx, Ex, gx = load(name)
    for s in range(3):
        Q = rand_codes(Vx, np.random.default_rng(100 + s))
        a0, b0 = leverages(Vx, nPx, Ex, Q)
        ld, err = flat_gauge(Vx, nPx, Ex, Q, 0, 3000, 1e-13)
        a, b = leverages(Vx, nPx, Ex, Q, ld)
        n0 = nullity(Vx, nPx, Ex, Q)[0]; n1 = nullity(Vx, nPx, Ex, Q, logd=ld)[0]
        check(err < 1e-10 and abs(np.sum(a * b) - Ex / 4) < 1e-8 and n0 == n1,
              f'C2 {name} seed{s}: raw lev range [{a0.min():.3f},{a0.max():.3f}] -> flat (err {err:.0e}); sum a_e b_e = {np.sum(a*b):.6f} = E/4; dim W {n0} invariant')
# K_{4,4}, girth 4: product code attains E/4
Vk = [[4 * i + j for j in range(4)] for i in range(4)] + [[4 * i + j for i in range(4)] for j in range(4)]
rng = np.random.default_rng(7); Ar = rand_codes([Vk[0]], rng)[0]; Bc = rand_codes([Vk[0]], rng)[0]
Qk = [Ar] * 4 + [Bc] * 4
nk = nullity(Vk, 4, 16, Qk)[0]
ldk, errk = flat_gauge(Vk, 4, 16, Qk, 0, 3000, 1e-13); ak, bk = leverages(Vk, 4, 16, Qk, ldk)
check(nk == 4 and abs(np.sum(ak * bk) - 4) < 1e-8, f'C2 K44 non-flat product code: dim W = {nk} = E/4, flat-gauge T1 = {np.sum(ak*bk):.6f} (tight)')
Qr = rand_codes(Vk, rng); check(nullity(Vk, 4, 16, Qr)[0] <= 4, 'C2 K44 random codes: dim W = %d <= 4' % nullity(Vk, 4, 16, Qr)[0])

# C3: bootstrap (field-free) bound and tightness on PG3
c, order = boot(V, E, R=300)
check(c == 8, f'C3 PG3 bootstrap bound dim W <= {c} ({c/E:.4f} < 1/6)')
f8 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sol_PG3_r8_s0.npy')
if os.path.exists(f8):
    Q8 = [orth(np.asarray(q, complex)) for q in np.load(f8, allow_pickle=True)]
    n8, w8 = nullity(V, nP, E, Q8, tol=1e-9)
    mg = min(abs(np.linalg.det(Y[list(S)])) for Y in Q8 for S in itertools.combinations(range(4), 2))
    check(n8 == 8 and w8[8] > 0.2 and mg > 0.1, f'C3 expXE r8 optimum: dim W = {n8}, spectral gap {w8[8]:.3f}, min 2x2 minor {mg:.3f} -> MDS, s = 8/52 = bootstrap (tight)')
for name in ['BA4', 'PG4M', 'W3', 'W4M', 'PG5']:
    Vx, nPx, Ex, gx = load(name); cx, _ = boot(Vx, Ex, R=100)
    check(cx / Ex <= 1 / gx + 1e-12, f'C3 {name}: bootstrap {cx}/{Ex} = {cx/Ex:.4f} <= 1/girth = {1/gx:.4f}  (closed in every field)')

# C4: gauge artefact
def perfect_matchings(Vx, nPx, Ex):
    from scipy.optimize import linear_sum_assignment
    ends = {}
    for v, es in enumerate(Vx):
        for e in es: ends.setdefault(e, []).append(v)
    left = list(range(nPx)); right = list(range(nPx, len(Vx))); avail = set(range(Ex)); Ms = []
    while avail:
        C = np.full((len(left), len(right)), 1e6)
        eid = {}
        for e in avail:
            p, l = ends[e]; C[p, l - nPx] = 0; eid[(p, l - nPx)] = e
        r, cidx = linear_sum_assignment(C); M = [eid[(i, j)] for i, j in zip(r, cidx)]
        assert C[r, cidx].sum() == 0; Ms.append(M); avail -= set(M)
    return Ms
Ms = perfect_matchings(V, nP, E); X = set(Ms[0]) | set(Ms[1])
Q = rand_codes(V, np.random.default_rng(11)); n0 = nullity(V, nP, E, Q)[0]
for t in [1e-1, 1e-2, 1e-3, 1e-4]:
    ld = np.array([0.0 if e in X else np.log(t) for e in range(E)])
    nl, w = nullity(V, nP, E, Q, tol=1e-30, logd=ld)
    print(f'     t={t:g}: eig_26(K) = {w[25]:.2e}, eig_27(K) = {w[26]:.2e}')
check(w[25] < 1e-5 and w[26] > 0.1 and n0 == 0, f'C4 random MDS code on PG3 has dim W = {n0}, yet 26 = E/2 eigenvalues of K -> 0 (like t^2) along the gauge orbit')

# C5: product
Vk_, Qk_ = Vk, Qk
PPa, PLa = projectors(Vk, 4, 16, Qk); PPb, PLb = projectors(V, nP, E, Q8 if os.path.exists(f8) else Qf)
K = 2 * np.eye(16 * E) - np.kron(PPa, PPb) - np.kron(PLa, PLb)
wk = np.linalg.eigvalsh(K); nprod = int(np.sum(wk < 1e-9))
ma = moments(PPa, PLa, 2); mb = moments(PPb, PLb, 2)
check(nprod == nk * (8 if os.path.exists(f8) else 0), f'C5 dim of product intersection = {nprod} = {nk} x dim W\'')
print('     product-certificate bound Theta(S)Theta(S\') =', round(markov_ref(ma, 2) * markov_ref(mb, 2), 6), '>= s(S)s(S\') =', round(nprod / (16 * E), 6))
print('PASS_EXPXG' if ok else 'FAILED')
