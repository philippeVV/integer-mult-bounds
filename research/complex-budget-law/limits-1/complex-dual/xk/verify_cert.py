"""Independent stress test of a saved local certificate h (split2 basis): many random-start maximisations of the
VERTEX residual and a fine grid for the EDGE residual."""
import numpy as np, sys
CERT = sys.argv[1]; sys.argv = [sys.argv[0], sys.argv[2], '6', '3']
src = open('split2.py').read().split('rng = np.random.default_rng')[0]
exec(src)
c = np.load(CERT)
rng = np.random.default_rng(12345)
for n in NS:
    best = -1
    for s in range(300):
        r = minimize(lambda z: -vval(proj(z, n), c), rng.standard_normal(n * n), method='BFGS', options=dict(maxiter=500))
        v = vval(proj(r.x, n), c)
        if v > best: best = v; a = np.sort(np.real(np.diag(proj(r.x, n))))
    print(f'n={n}: max VERTEX residual per coordinate over 300 starts = {best:.3e} (target scale n/64 per vertex = {1/64:.4f}/coord) at lev {np.round(a,3)}', flush=True)
xs = np.linspace(0, 1, 401); mx = (-1, None)
for x in xs:
    m = x * (1 - x); us = np.linspace(0, m, 81); U, Wg = np.meshgrid(us, us)
    viol = 3 * U * Wg - e(x) - feat(x, U) @ c - feat(1 - x, Wg) @ c
    if viol.max() > mx[0]: mx = (viol.max(), x)
print('EDGE residual max on 401x81x81 grid:', mx)
