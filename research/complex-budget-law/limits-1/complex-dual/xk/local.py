"""Local (tree) formulas for m1, m2, m3 at girth >= 8 in terms of the per-vertex projectors and the opposite leverages."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from xk_sparse import *

def local_moments(S, Qs, x=None):
    a, b = leverages(S, Qs, x); E = S.E
    m1 = np.sum(a * b); m2 = m3 = 0.0; u = np.zeros(E); w = np.zeros(E)
    for i, (es, Q) in enumerate(zip(S.V, Qs)):
        Y = orth(Q if x is None else np.exp(x[es])[:, None] * Q); P = Y @ Y.conj().T
        d = (b if i < S.nP else a)[es]; PD = P * d[None, :]
        m2 += np.real(np.trace(PD @ PD)); m3 += np.real(np.trace(PD @ PD @ PD))
        off = abs(P) ** 2; np.fill_diagonal(off, 0)
        (u if i < S.nP else w)[es] = off @ d
    m2 -= np.sum(a**2 * b**2); m3 += -np.sum(a**3 * b**3) + 3 * np.sum(u * w)
    return np.array([1, m1 / E, m2 / E, m3 / E]), a, b, u, w

if __name__ == '__main__':
    for name in ['W3', 'W5']:
        V, nP, E, g = load_any(name); S = Struct(V, nP, E); rng = np.random.default_rng(1)
        Q = rand_codes(V, rng); x = rng.standard_normal(E) * 0.7
        Ys = [orth(np.exp(x[es])[:, None] * q) for es, q in zip(V, Q)]
        print(name, theta_grad(S, Ys, 3)[1], local_moments(S, Q, x)[0])
