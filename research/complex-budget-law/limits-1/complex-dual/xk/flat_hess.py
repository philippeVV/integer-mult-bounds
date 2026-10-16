"""Second-order behaviour of the vertex functional Psi_0 at flat-leverage projectors (eta drops out at 2nd order)."""
import numpy as np, sys
sys.argv = sys.argv[:1] + ['0.5', '4', '7', sys.argv[1] if len(sys.argv) > 1 else '0']
exec(open('split_lp.py').read().split('pool = []')[0])
for n in [4, 6]:
    k = n // 2; F = np.exp(2j*np.pi*np.outer(np.arange(n), np.arange(n))/n)/np.sqrt(n)
    rng = np.random.default_rng(1); worst = -1
    for trial in range(200):
        Y0 = F[:, rng.permutation(n)[:k]] if trial % 2 else np.linalg.qr(rng.standard_normal((n,k)) + 1j*rng.standard_normal((n,k)))[0]
        P0 = Y0 @ Y0.conj().T
        if abs(np.real(np.diag(P0)) - 0.5).max() > 1e-9: continue
        D = rng.standard_normal((n, k)) + 1j*rng.standard_normal((n, k))
        for h in [1e-2, 3e-3]:
            Y = np.linalg.qr(Y0 + h * D)[0]; P = Y @ Y.conj().T
            v = parts(P)[0] / h**2; worst = max(worst, v)
    print('n', n, 'max second-order coeff of Psi_0 - target over random directions at flat points:', worst)
