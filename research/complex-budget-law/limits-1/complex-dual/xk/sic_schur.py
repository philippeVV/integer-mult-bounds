"""Spectrum of the code/gauge Schur complement at SIC codes on W3: count zero modes vs the flat-manifold tangent dimension."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from maxmin_hess import *
from sic import sic_codes
V, nP, E, g = load_any('W3'); S = Struct(V, nP, E); rng = np.random.default_rng(0)
Qs = sic_codes(V, rng)
f0, g0, H, nw = hessian(S, Qs)
Hww, Hwy, Hyy = H[:nw, :nw], H[:nw, nw:], H[nw:, nw:]
Sm = Hww - Hwy @ np.linalg.pinv(Hyy, rcond=1e-8) @ Hwy.T; ev = np.linalg.eigvalsh((Sm + Sm.T) / 2)
sc = abs(ev).max()
neg = ev[ev < -1e-7 * abs(ev).max()]; print('largest of the strictly negative Schur eigenvalues: %.3e' % neg.max()); print('code dims', nw, ' Schur eig: max %.2e, #>-1e-7*scale: %d, #<-1e-7*scale: %d, scale %.2e' % (ev[-1], np.sum(ev > -1e-7 * sc), np.sum(ev < -1e-7 * sc), sc))
# flat-manifold tangent: code directions with zero first-order leverage change on both sides
perp = setup(V, Qs); rows = []
for i, (es, Y, Pp) in enumerate(zip(V, Qs, perp)):
    n, k = Y.shape; m = (n - k) * k
    for j in range(2 * m):
        wv = np.zeros(2 * m); wv[j] = 1; W = (wv[:m] + 1j * wv[m:]).reshape(n - k, k)
        dY = Pp @ W; dlev = 2 * np.real(np.sum(dY * Y.conj(), axis=1))
        rows.append((i, j, dlev))
offs = np.cumsum([0] + [2 * (Y.shape[0] - Y.shape[1]) * Y.shape[1] for Y in Qs])
Jm = np.zeros((len(V) * 4, nw))
for (i, j, dl) in rows: Jm[4 * i: 4 * i + 4, offs[i] + j] = dl
r = np.linalg.matrix_rank(Jm, 1e-9); print('flat-manifold tangent dim =', nw - r, ' (code dims', nw, ', leverage-constraint rank', r, ')')
# Schur restricted to the complement of the flat tangent
U, s, Vt = np.linalg.svd(Jm); Cmp = Vt[:r].T   # row space of Jm = complement of its kernel
ev2 = np.linalg.eigvalsh(Cmp.T @ Sm @ Cmp); print('Schur on complement of flat tangent: max eig %.3e, min %.3e' % (ev2[-1], ev2[0]))
# packaging check (limits-1): 559 zero modes = 400 flat-tangent + 159 gauge; the 81 normal directions are strictly negative
assert np.sum(ev > -1e-7 * sc) == 559 and np.sum(ev < -1e-7 * sc) == 81 and nw - r == 400
assert ev2[-1] < 1e-8 and np.sum(ev2 < -1e-7 * sc) == 81
print('PASS_EXPXK_SCHUR')
