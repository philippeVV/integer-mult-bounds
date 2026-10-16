"""Global test along rays delta -> t*delta at flat degree-4 codes: linear (Schur) gauge rule vs true inner inf."""
from model import *
from xk_sparse import inner_min
name = sys.argv[1]; V, nP, E, g = load_any(name); S = Struct(V, nP, E); rng = np.random.default_rng(int(sys.argv[2]))
conc = float(sys.argv[3])
W = rng.dirichlet([conc]*3, size=len(V)) if conc > 0 else np.full((len(V), 3), 1/3)
Q0 = [flat_code(W[i], rng) for i in range(len(V))]
X, Y = XYmats(V, nP, E, W); K = Kmat(X, Y); G, C, D = xdelta(K); ev, es_, Sm = schur(K)
print('Gamma min nontriv', np.sort(ev)[1], 'Schur max', es_[-1])
Gp = np.linalg.pinv(G, rcond=1e-10)
for trial in range(4):
    if trial < 2: d = rng.standard_normal(E)
    else:  # least negative Schur direction (excluding trivial vertex-constant kernel)
        w_, U = np.linalg.eigh(Sm); nz = np.where(w_ < -1e-9 * abs(w_).max())[0]; d = U[:, nz[-1 - (trial - 2) * 5]]
    d /= np.linalg.norm(d) / np.sqrt(E)
    xs = -Gp @ C @ d
    out = []
    for t in [0.05, 0.2, 0.5, 1, 2, 4]:
        Qs = [Q0[i] * (np.exp(t * d[es])[:, None] if i >= nP else 1) for i, es in enumerate(V)]
        lin = theta(S, Qs, 3, t * xs)
        inf, _ = inner_min(S, Qs, 3, [t * xs, np.zeros(E)])
        out.append(f't={t}: lin {lin:.6f} inf {inf:.6f} model {0.125 + 4/E*t*t*(d@Sm@d):.6f}')
    print('dir', trial, ' | '.join(out), flush=True)
