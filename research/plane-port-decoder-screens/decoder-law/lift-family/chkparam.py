import numpy as np, sys
from param import *
D = np.load('sols_0_20.npz')
for key in sorted(D.keys(), key=lambda k: tuple(map(int, k.split('_')))):
    V, W = D[key]; N = V @ W.T
    dev = odd & (abs(abs(N * N.T) - 0.25) > 1e-6); deg = dev.sum(1)
    S = [a for a in range(NL) if deg[a] >= 50]
    T = [a for a in range(NL) if a not in S]
    out = [key, 'S', len(S), 'pairwise-odd', all(odd[a, b] for a in S for b in S if a != b)]
    for side, (X, Y) in (('X', (V, W)), ('Y', (W, V))):
        rT = np.linalg.matrix_rank(X[T], 1e-7)
        if rT != 8: out += [side, 'rankT', rT]; continue
        X2, Y2, d, err, fit = gauge_fix(X, Y, S)
        Ss = np.array(S, int); Ta = np.array(T)
        xi = X2[Ss, 8]; a = X2[Ss, :8]
        c = a / xi[:, None] if len(S) else np.zeros((0, 8))
        # Y: y_t[:8] = r_t/2 ?  y_s[:8] = beta_s r_s ?
        e1 = np.abs(Y2[Ta, :8] - R[Ta] / 2).max()
        beta = np.array([Y2[s, :8] @ R[s] / 2 for s in S]); e2 = max([np.abs(Y2[s, :8] - beta[i] * R[s]).max() for i, s in enumerate(S)] + [0])
        eta = Y2[:, 8]
        e3 = max([abs(eta[t] + c[i] @ R[t] / 2) for i, s in enumerate(S) for t in T if IP[s][t] == 0] + [0])
        # is eta linear on T?  residual of lstsq eta_T = R_T lam
        lam = np.linalg.lstsq(R[Ta], eta[Ta], rcond=None)[0]; e4 = np.abs(R[Ta] @ lam - eta[Ta]).max()
        F = [t for t in T if all(IP[s][t] != 0 for s in S)]
        # c_s = c0 + mu_s r_s ?
        if len(S) >= 1:
            # fit c0 common: c_s - c0 in span r_s ;  solve least squares for c0 and mu
            nS = len(S); Am = np.zeros((8 * nS, 8 + nS)); b = c.reshape(-1)
            for i, s in enumerate(S): Am[8*i:8*i+8, :8] = np.eye(8); Am[8*i:8*i+8, 8+i] = R[s]
            sol = np.linalg.lstsq(Am, b, rcond=None)[0]; e5 = np.abs(Am @ sol - b).max()
        else: e5 = 0
        out += [side, f'torusErr {err:.1e} fit {fit:.1e} yT {e1:.1e} ys {e2:.1e} etaOrth {e3:.1e} etaLin {e4:.1e} c=c0+mu r {e5:.1e} |F| {len(F)} minxi {np.abs(xi).min() if len(S) else 0:.2e} minbeta {np.abs(beta).min() if len(S) else 0:.2e}']
    print(*out, flush=True)
