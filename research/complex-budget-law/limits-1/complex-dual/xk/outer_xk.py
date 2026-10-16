"""J=3 adversary (sparse): maximise over local codes  inf_{|log d|<=BOX} Theta_3(D codes).
usage: outer_xk.py GRAPH seed steps init(rand|fourier|file.npy) BOX [eta] [maxiter]"""
import sys, os, json, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from xk_sparse import *
name, seed, steps, init, BOX = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], float(sys.argv[5])
eta = float(sys.argv[6]) if len(sys.argv) > 6 else 0.03
MI = int(sys.argv[7]) if len(sys.argv) > 7 else 200
J = 3
V, nP, E, g = load_any(name); S = Struct(V, nP, E); rng = np.random.default_rng(seed)
if init == 'rand': Qs = rand_codes(V, rng)
elif init == 'fourier': Qs = [orth(Q + 0.05 * (rng.standard_normal(Q.shape) + 1j * rng.standard_normal(Q.shape))) for Q in fourier_codes(V, rng)]
else: Qs = [orth(np.asarray(Q, complex)) for Q in np.load(init, allow_pickle=True)]
tag = f'{name}_box{BOX:g}_s{seed}_{os.path.basename(init).split(".")[0]}'
x = np.zeros(E); best = (-1, None, None); t0 = time.time()
for it in range(steps):
    starts = [x, np.zeros(E), flat_gauge(S, Qs, 0, 300)[0], flat_gauge(S, Qs, 1, 300)[0]]
    if it % 5 == 0: starts.append(rng.standard_normal(E))
    val, x = inner_min(S, Qs, J, starts, box=BOX, maxiter=MI)
    interior = float(abs(x).max()) < BOX - 1e-3
    if val > best[0] and interior: best = (val, [Q.copy() for Q in Qs], x.copy())
    Ys = [orth(np.exp(x[es])[:, None] * Q) for es, Q in zip(V, Qs)]
    th, Zs, m = code_grad(S, Ys, J)
    gn = np.sqrt(sum(np.sum(abs(Z) ** 2) for Z in Zs)) + 1e-15
    Qs = [orth(np.exp(-x[es])[:, None] * (Y + eta * Z / gn)) for es, Y, Z in zip(V, Ys, Zs)]
    a, b = leverages(S, Ys)
    print(json.dumps(dict(it=it, theta_star=round(val, 5), m=[round(v, 5) for v in m[1:]], gn=round(gn, 5),
          levP=[round(a.min(), 3), round(a.max(), 3)], levL=[round(b.min(), 3), round(b.max(), 3)], xmax=round(float(abs(x).max()), 2),
          interior=interior, best=round(best[0], 5), t=round(time.time() - t0))), flush=True)
    if best[1] is not None and it % 10 == 9:
        np.save(f'adv_{tag}.npy', np.array(best[1], dtype=object), allow_pickle=True); np.save(f'adv_{tag}_gauge.npy', best[2])
val, Qb, xb = best
print(json.dumps(dict(graph=name, J=J, box=BOX, best_interior_theta_star=val, threshold=1 / g)))
if Qb is not None:
    np.save(f'adv_{tag}.npy', np.array(Qb, dtype=object), allow_pickle=True); np.save(f'adv_{tag}_gauge.npy', xb)
np.save(f'adv_{tag}_last.npy', np.array(Qs, dtype=object), allow_pickle=True)
