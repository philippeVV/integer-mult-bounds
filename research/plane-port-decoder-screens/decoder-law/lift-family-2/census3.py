"""All 28 saved completions (expSJ/sols_0_20.npz): find the hidden lift set A (support of the sparsest X-deficient line,
which has rank-8 complement), fit the kappa-family normal form (expSM sec 0) with S = A, read B = supp kappa', and compare
the moduli tangent dimension (tangent.py) with the kappa-family dimension |A| + |B| - 1 - [|A| = 1]."""
import sys, numpy as np, itertools
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from e8 import *
from param import gauge_fix
from tangent import jac, orbit, nrank
D = np.load('./sols_0_20.npz'); R = np.array(ROOT, float) / 2
CLS = {h: np.array([dotq(LAB[a], h) for a in range(NL)], bool) for h in range(1, 511)}
def sparsest(M):
    best = None
    for h, O in CLS.items():
        u, s, vt = np.linalg.svd(M[O], full_matrices=False)
        if s[-1] < 1e-7 * s[0]:
            v = M @ vt[-1]; sup = list(np.where(np.abs(v) > 1e-7 * np.abs(v).max())[0])
            if best is None or len(sup) < len(best): best = sup
    return best
allok = True
for k in sorted(D.keys(), key=lambda k: tuple(map(int, k.split('_')))):
    V, W = D[k]
    out = None
    cands = []
    for side, (X, Y) in (('X', (V, W)), ('Y', (W, V))):
        A = sparsest(X); T = [a for a in range(NL) if a not in A]
        if A is None or np.linalg.matrix_rank(X[T], 1e-7) != 8: continue
        cands.append((len(A), side, X, Y, A))
    for _, side, X, Y, A in sorted(cands, key=lambda c: c[0])[:1]:
        X2, Y2, d, err, fit = gauge_fix(X, Y, A)
        # c0 gauge (expSM sec 0 (ii)): shift so that x_s[:8] - ka_s c0 is parallel to r_s for all s in A
        ka0 = X2[A, 8]; rows = []; rhs = []
        for i, s_ in enumerate(A):
            Pp = np.eye(8) - np.outer(R[s_], R[s_]) / 2
            rows.append(ka0[i] * Pp); rhs.append(Pp @ X2[s_, :8])
        c0 = np.linalg.lstsq(np.vstack(rows), np.concatenate(rhs), rcond=None)[0]
        Mg = np.eye(9); Mg[8, :8] = -c0
        X2 = X2 @ Mg; Y2 = Y2 @ np.linalg.inv(Mg).T
        al = np.array([X2[a, :8] @ R[a] / 2 for a in range(NL)]); alp = np.array([Y2[a, :8] @ R[a] / 2 for a in range(NL)])
        res = max(max(np.abs(X2[a, :8] - al[a] * R[a]).max(), np.abs(Y2[a, :8] - alp[a] * R[a]).max()) for a in range(NL))
        ka, kap = X2[:, 8], Y2[:, 8]
        diag = np.abs(2 * al * alp + ka * kap - 1).max()
        sc = np.abs(kap).max()
        B = list(np.where(np.abs(kap) > 1e-7 * sc)[0])
        suppA_ok = set(np.where(np.abs(ka) > 1e-7 * np.abs(ka).max())[0]) == set(A)
        pw = all(IP[a][b] != 0 for a, b in itertools.combinations(A, 2)); cross = all(IP[a][b] != 0 for a in A for b in B if a != b)
        out = (side, A, B, res, diag, suppA_ok, pw, cross); break
    J = jac(V, W); rJ, _ = nrank(J, ''); O = orbit(V, W); rO, _ = nrank(O, '')
    mod = 2 * NL * 9 - rJ - rO
    if out is None:
        print(k, 'NO LIFT FIT', 'moduli', mod); allok = False; continue
    side, A, B, res, diag, sok, pw, cross = out
    pred = len(A) + len(B) - 1 - (len(A) == 1)
    ok = res < 1e-9 and diag < 1e-9 and sok and pw and cross and pred == mod
    a2 = len(A) == 3 and np.linalg.matrix_rank(R[A]) == 2
    allok &= ok
    print(k, 'rigid side', side, '|A|', len(A), '|B|', len(B), 'A=B' if set(A) == set(B) else '', f'normal-form residual {res:.1e} diag {diag:.1e}',
          'A pairwise odd', pw, 'A-B non-orth', cross, 'A2-triple' if a2 else '', 'pred dim', pred, 'moduli tangent', mod, 'OK' if ok else 'CHECK', flush=True)
print('ALL LIFTS, tangent = kappa-family dim' if allok else 'SOME CHECK')
