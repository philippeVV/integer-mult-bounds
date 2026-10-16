"""Stop-test adversary for Conjecture M: maximise J(d) = Phi(d) - Phi(lam d) on the sphere |d|_rms = r (J > 0 <=> non-monotone ray).
Danskin: grad Phi(d) = g_v at the optimal gauge.  Also reports the exact radial derivative R = <g_u,u> + <g_v,v> = Phi'(1) at t = 1.
usage: m_adv.py GRAPH base seed r steps [lam]"""
import sys, time
from xr_tools import *
from xk_sparse import flat_gauge
name, base, seed, r, steps = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
lam = float(sys.argv[6]) if len(sys.argv) > 6 else 0.8
S, Q0, W, G, C, D, lmin, rng = setup(name, base, seed); E = S.E
print(f'{name} base {base} seed {seed} r {r} lam {lam}: lambda_min(Gamma) {lmin:.3e}', flush=True)
Lr = -np.linalg.pinv(G, rcond=1e-10) @ C
Bv = np.zeros((E, len(S.V)))
for i, es in enumerate(S.V): Bv[es, i] = 1
Ub, sb, _ = np.linalg.svd(Bv, full_matrices=False); Ub = Ub[:, sb > 1e-9]     # span of vertex constants (trivial delta)
def proj(d): return d - Ub @ (Ub.T @ d)
def norm(d): d = proj(d); return d * r * np.sqrt(E) / np.linalg.norm(d)
xs = {}
def Phi(d, key):
    other = xs.get(1 - key)
    st = ([xs[key]] if key in xs else []) + [Lr @ d, -d / 2, np.zeros(E)]
    if other is not None: st.append(other * (lam if key == 0 else 1 / lam))
    th, x, gu, gv = inner(S, Q0, d, st); xs[key] = x
    return th, x, gu, gv
d = norm(rng.standard_normal(E)); eta = 0.3; t0 = time.time()
def J(d):
    a = Phi(d, 1); b = Phi(lam * d, 0)
    return a[0] - b[0], a[3] - lam * b[3], a, b
j, gJ, a, b = J(d); best = (j, d.copy())
for it in range(steps):
    gq = proj(gJ); gt = gq - (gq @ d) / (d @ d) * d          # tangent to the sphere
    if np.linalg.norm(gt) < 1e-14: break
    dn = norm(d + eta * r * np.sqrt(E) * gt / np.linalg.norm(gt))
    jn, gJn, an, bn = J(dn)
    if jn > j: d, j, gJ, a, b = dn, jn, gJn, an, bn; eta = min(eta * 1.3, 1.0)
    else: eta *= 0.5
    if j > best[0]: best = (j, d.copy())
    x = a[1]; R = a[2] @ x + a[3] @ (x + d)
    if it % 5 == 0 or it == steps - 1:
        print(f'it {it}: J = {j:.3e}  Phi(1) = {a[0]:.9f}  Phi(lam) = {b[0]:.9f}  R = Phi\'(1) = {R:.3e}  |grad_x| = {np.abs(a[2]+a[3]).max():.1e}  eta {eta:.3g}  t {time.time()-t0:.0f}', flush=True)
    if eta < 1e-6: break
np.save(f'madv_{name}_{base}_{seed}_r{r:g}.npy', best[1])
# precise re-solve of the best ray: forward and backward continuation + fresh starts, min over all
d = best[1]; ts = np.linspace(0.1, 1.5, 15); vals = np.full(len(ts), np.inf); xw = {}
for sweep in [range(len(ts)), range(len(ts) - 1, -1, -1)]:
    xp = None
    for i in sweep:
        t = ts[i]; st = [Lr @ (t * d), -t * d / 2, np.zeros(E)] + ([xp] if xp is not None else [])
        th, xp, gu, gv = inner(S, Q0, t * d, st); vals[i] = min(vals[i], th)
mono = all(vals[i + 1] <= vals[i] + 1e-10 for i in range(len(ts) - 1))
print('ray re-solve Phi(t), t=0.1..1.5:', ' '.join(f'{v:.9f}' for v in vals), 'monotone' if mono else 'NON-MONOTONE', f'max {vals.max():.9f}')
print(f'SUMMARY M {name} base {base} r {r}: max J = {best[0]:.3e} ({"NON-MONOTONE candidate" if best[0] > 1e-9 else "monotone"})')
