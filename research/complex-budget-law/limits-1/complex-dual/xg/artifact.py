"""Gauge-artifact check on expXE saved optima (read-only): their loss (r smallest eigenvalues of K=2I-Pi_P-Pi_L)
is not gauge invariant; dim W is.  Re-balance the gauge (Barthe flat on P, then on L) and recount."""
import sys, numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from xg_core import *
V, nP, E, g = load('PG3')
for f in sys.argv[1:]:
    Qs = list(np.load(f, allow_pickle=True))
    Qs = [orth(np.asarray(Q, complex)) for Q in Qs]
    n0, w0 = nullity(V, nP, E, Qs)
    a, b = leverages(V, nP, E, Qs)
    ld, err = flat_gauge(V, nP, E, Qs, 0, 5000, 1e-13)
    n1, w1 = nullity(V, nP, E, Qs, logd=ld)
    a1, b1 = leverages(V, nP, E, Qs, ld)
    # exact-ish nullity: singular values of the stacked constraint matrix [A_P^perp ; A_L^perp] (gauge-free count)
    rows = []
    for es, Q in zip(V, Qs):
        n = len(es); Pperp = np.eye(n) - Q @ Q.conj().T
        U, sv, _ = np.linalg.svd(Pperp); C = np.zeros((n - n // 2, E), complex); C[:, es] = U[:, : n - n // 2].conj().T
        rows.append(C)
    sv = np.linalg.svd(np.vstack(rows), compute_uv=False)
    print(f, '\n raw: 12 smallest eig(K)', np.round(w0[:12], 9), ' lev P range', a.min().round(4), a.max().round(4),
          '\n P-flat gauge (err %.1e): 12 smallest eig(K)' % err, np.round(w1[:12], 6), ' lev L range', b1.min().round(4), b1.max().round(4),
          '\n constraint-matrix smallest singular values', np.round(sv[-12:], 9), 'cond-free nullity(1e-6 rel)', int(np.sum(sv < 1e-6 * sv[0])))
