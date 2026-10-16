"""Initial codes: SIC (degree 4) / icosahedral ETF(3,6) (degree 6), random relabelling + phases + small noise."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from sic import SIC
from xk_sparse import load_any, orth
name, noise, seed = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]); V, nP, E, g = load_any(name); rng = np.random.default_rng(seed)
phi = (1 + 5 ** 0.5) / 2
ico = np.array([[0, 1, phi], [0, 1, -phi], [1, phi, 0], [1, -phi, 0], [phi, 0, 1], [-phi, 0, 1]], float)
ico /= np.linalg.norm(ico, axis=1, keepdims=True); ICO = ico / np.sqrt(2)   # 6 unit vectors in R^3, frame op = 2 I
out = []
for es in V:
    base = SIC if len(es) == 4 else ICO
    Y = base[rng.permutation(len(es))] * np.exp(2j * np.pi * rng.random(len(es)))[:, None]
    Y = Y + noise * (rng.standard_normal(Y.shape) + 1j * rng.standard_normal(Y.shape))
    out.append(orth(Y))
P = out[0] @ out[0].conj().T; print('|P|^2 offdiag range', np.round(np.sort(np.unique(np.round(abs(P[np.triu_indices(len(P),1)])**2, 4))), 4), 'diag', np.round(np.real(np.diag(P)), 4))
np.save(f'init_{name}_etf_n{noise:g}_s{seed}.npy', np.array(out, dtype=object), allow_pickle=True)
