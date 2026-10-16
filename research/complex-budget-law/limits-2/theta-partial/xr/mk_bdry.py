"""Starts on the boundary of {Gamma >= 0}: bisect w between SIC (Gamma > 0) and a random w (Gamma indefinite);
add a small relative gauge delta along the least-negative Schur direction. Saves codes for ../expXK/outer_xk.py."""
from model import *
name, seed, conc, eps = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
V, nP, E, g = load_any(name); rng = np.random.default_rng(seed)
W0 = np.full((len(V), 3), 1/3); W1 = rng.dirichlet([conc]*3, size=len(V))
def lam(W):
    X, Y = XYmats(V, nP, E, W); K = Kmat(X, Y); G, C, D = xdelta(K); o = np.ones(E)/np.sqrt(E); P = np.eye(E)-np.outer(o,o); return np.linalg.eigvalsh(P@G@P + np.outer(o,o))[0], K
l0, _ = lam(W0); l1, _ = lam(W1); print('lam SIC', l0, 'lam rand', l1)
if l0 <= 0: print('SIC not Gamma>0 on this graph'); sys.exit()
if l1 > 0: print('random w already Gamma>0; using it'); lo = 1.0
else:
    lo, hi = 0.0, 1.0
    for _ in range(40):
        mid = (lo + hi) / 2
        if lam((1 - mid) * W0 + mid * W1)[0] > 0: lo = mid
        else: hi = mid
W = (1 - lo) * W0 + lo * W1; l, K = lam(W); ev, es_, Sm = schur(K)
w_, U = np.linalg.eigh(Sm); nz = np.where(w_ < -1e-9 * abs(w_).max())[0]; d = U[:, nz[-1]]; d /= np.linalg.norm(d) / np.sqrt(E)
print('boundary s', lo, 'lam_min', l, 'least-negative Schur eig', w_[nz[-1]])
Qs = [flat_code(W[i], rng) * (np.exp(eps * d[es])[:, None] if i >= nP else 1) * np.exp(2j*np.pi*rng.random(4))[:, None] for i, es in enumerate(V)]
Qs = [orth(q) for q in Qs]
np.save(f'bdry_{name}_s{seed}_e{eps:g}.npy', np.array(Qs, dtype=object), allow_pickle=True)
