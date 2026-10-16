"""Ray monotonicity test: Phi(t) = inf_x Theta3(flat(w) points scaled e^x, lines scaled e^{x + t delta}).
usage: ray2.py GRAPH seed conc ndirs"""
from model import *
from check_xr_util import XY_codes
from xk_sparse import inner_min, flat_gauge
from comp import comp_gauge
name, seed, conc, nd = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
V, nP, E, g = load_any(name); S = Struct(V, nP, E); rng = np.random.default_rng(seed)
W = rng.dirichlet([conc]*3, size=len(V)) if conc > 0 else np.full((len(V), 3), 1/3)
from check_xr_util import XY_codes
if conc >= 0:
    Q0 = [flat_code(W[i], rng) for i in range(len(V))]
else:   # conc < 0: generic flat codes (any degree): Barthe-flattened random complex codes
    from check_xr_util import flatten
    import xg_core as Dc
    Q0 = [flatten(q) for q in Dc.rand_codes(V, rng)]
X, Y = XY_codes(V, nP, E, Q0); K = Kmat(X, Y); G, C, D = xdelta(K)
o = np.ones(E)/np.sqrt(E); lmin = np.linalg.eigvalsh((np.eye(E)-np.outer(o,o))@G@(np.eye(E)-np.outer(o,o)) + np.outer(o,o))[0]
print(f'{name} seed {seed} conc {conc}: lambda_min(Gamma) = {lmin:.3e}', flush=True)
ts = [0, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6] if E < 500 else [0, 0.1, 0.3, 0.6, 1, 1.5, 2.5, 4]
gmax = 0; nonmono = 0
for k in range(nd):
    d = rng.standard_normal(E) if k % 2 == 0 else rng.standard_normal(E) * (rng.random(E) < 0.15) * 3
    d /= np.linalg.norm(d) / np.sqrt(E); x = np.zeros(E); vals = []
    for t in ts:
        Qs = [Q0[i] * (np.exp(t * d[es])[:, None] if i >= nP else 1) for i, es in enumerate(V)]
        st = [x, np.zeros(E), -t * d / 2, flat_gauge(S, Qs, 0, 300)[0], flat_gauge(S, Qs, 1, 300)[0]]
        try:
            xc = comp_gauge(S, Qs, x)[0]
            if np.all(np.isfinite(xc)): st.append(xc)
        except Exception: pass
        val, x = inner_min(S, Qs, 3, st, maxiter=3000); vals.append(val)
    mono = all(vals[i+1] <= vals[i] + 1e-7 for i in range(len(vals)-1)); nonmono += not mono; gmax = max(gmax, max(vals))
    print(f'dir {k}: ' + ' '.join(f'{v:.5f}' for v in vals) + ('' if mono else '  NONMONO'), flush=True)
print(f'SUMMARY {name} seed {seed} conc {conc}: max Phi = {gmax:.7f}, non-monotone rays {nonmono}/{nd}')
