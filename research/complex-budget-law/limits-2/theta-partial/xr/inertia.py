from model import *
name = sys.argv[1]; V, nP, E, g = load_any(name); rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv)>2 else 0)
def inert(M, tol=1e-9):
    ev = np.linalg.eigvalsh(M); s = abs(ev).max(); return (int((ev > tol*s).sum()), int((ev < -tol*s).sum()), int((abs(ev) <= tol*s).sum()))
for conc in [0.3, 1, 3, 30]:
    for t in range(2):
        W = rng.dirichlet([conc]*3, size=len(V)); X, Y = XYmats(V, nP, E, W); K = Kmat(X, Y)
        G, C, D = xdelta(K); ev, es_, Sm = schur(K)
        print(f'conc {conc}: K inertia {inert(K)}  Gamma inertia {inert(G)}  Schur max {es_[-1]:.2e} min {es_[0]:.2e}  Gamma min {np.sort(ev)[1]:.3e}')
