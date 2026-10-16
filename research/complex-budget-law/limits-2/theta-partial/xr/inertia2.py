from model import *
name = sys.argv[1]; V, nP, E, g = load_any(name); rng = np.random.default_rng(1)
def inert(M, tol=1e-9):
    ev = np.linalg.eigvalsh((M+M.T)/2); s = abs(ev).max(); return (int((ev > tol*s).sum()), int((ev < -tol*s).sum()), int((abs(ev) <= tol*s).sum()))
I = np.eye(E)
for conc in [0.3, 3]:
    W = rng.dirichlet([conc]*3, size=len(V)); X, Y = XYmats(V, nP, E, W)
    H0 = np.block([[4*Y, 3*(I+Y@X)], [3*(I+X@Y), 4*X]])
    print('H0 full inertia', inert(H0))
    A1 = np.block([[I, I], [Y, X]]); A2 = np.block([[I, -I], [Y, -X]])
    print('rank sigma-phi', np.linalg.matrix_rank(A1), 'rank d-psi', np.linalg.matrix_rank(A2))
    # restrict to R_P + R_L
    Z = np.block([[I-X, 0*I], [0*I, I-Y]]); U, s, _ = np.linalg.svd(Z); R = U[:, s > 1e-9]
    print('restricted H0 inertia', inert(R.T@H0@R), 'rank sigma-phi on R', np.linalg.matrix_rank(A1@R), 'rank d-psi on R', np.linalg.matrix_rank(A2@R))
    ex, ey = np.linalg.eigvalsh(X), np.linalg.eigvalsh(Y); print('X spec range', ex[[0,-1]], 'mult of 1', (abs(ex-1)<1e-9).sum())
