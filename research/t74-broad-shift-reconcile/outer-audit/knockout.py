"""Per-layer kappa value by knockout: sup kappa over (eps, beta, xi) of the PR2044 outer row system with one layer's
row reverted to its predecessor form. Slacks r0, zeta, layout slack set to 0 (sup = boundary of the strict system).
Active set from a float LP (scipy), then the vertex is solved and every constraint re-checked in exact Fractions at the
package's full-precision q. Written from the row tables in FINE-PRIME-PADDING.md sec.7, UNION-SCHEDULE.md, GAUGEFREE-INTERIORS.md
sec.6 and PR1673 check_frontier_fusion.py; imports nothing from the PR."""
import json, sys, itertools
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
sys.set_int_max_str_digits(0)
P = sys.argv[1].rstrip('/') + '/'
spec = json.load(open(P + 'new/SPEC.json'))
DIG = spec['outer_refinement']['endpoint_digits']
u = F(1, 10**(DIG + 10))
cl = F(spec['expected']['complex_lower']); b = cl - u; q = b - u; a = (q + b) / 2
online = F(spec['expected']['kappa'])
# variables (eps, beta, xi, kappa); constraint form: coef . v + const >= 0
def build(cfg):
    c = cfg['c']; H, C = cfg['halo'], cfg['core']
    K = []
    def ge(name, e, bb, x, k, const): K.append((name, (F(e), F(bb), F(x), F(k)), F(const)))
    ge('main eps q', q, 0, 0, -1, 0)
    ge('router a', 0, 0, 0, -1, a)
    ge('gaussian 1-eps+q beta', -1, q, 0, -1, 1)
    ge('group bit router', -1, a, 0, -1, 1)
    ge('group scans', -1, 1, 0, -1, 1)
    ge('boundary c+xi-eps', -1, 0, 1, -1, c)
    if cfg['union'] == 'radix':
        ge('union xi+1/2-eps', -1, 0, 1, -1, F(1, 2))
    else:  # PR1673 per-axis-scan union row 1-eps-beta-1/2+xi
        ge('old union 1-eps-beta-1/2+xi', -1, -1, 1, -1, F(1, 2))
        if cfg.get('full_key'):
            ge('old union minus key', 0, -2, 1, -1, F(-1, 2))
    ge('polynomial eps', 1, 0, 0, -1, 0)
    ge('microbox face', -1, -1, 0, -1, 1 - (H - C + 2))
    ge('eps<=1', -1, 0, 0, 0, 1)
    ge('xi>=1/2', 0, 0, 1, 0, F(-1, 2))
    ge('xi<=eps', 1, 0, -1, 0, 0)
    ge('beta>=0', 0, 1, 0, 0, 0)
    ge('beta<=eps', 1, -1, 0, 0, 0)
    ge('beta<=1', 0, -1, 0, 0, 1)
    g = cfg['guard']
    if g == 'fine':      ge('guard 1-(2xi-2)-beta', 0, -1, -2, 0, 3)
    elif g == 'theta_d': ge('guard 1-(2xi-eps)-beta', 1, -1, -2, 0, 1)
    elif g == 'linear':  ge('guard 1-xi-beta', 0, -1, -1, 0, 1)
    ge('membership 1-beta-2(1-eps)', 2, -1, 0, 0, -1)
    ge('halo volume', 0, -1, 0, 0, C - H)
    ge('short record eps>=a', 1, 0, 0, 0, -a)
    if cfg.get('tag2'): ge('old tag 1-2beta', 0, -2, 0, 0, 1)
    if cfg.get('old_scan'): ge('old per-axis scan 1-eps', -1, 0, 0, -1, 1)
    return K
def solve(cfg):
    K = build(cfg)
    A = np.array([[-float(x) for x in co] for _, co, _ in K]); bb = np.array([float(k) for _, _, k in K])
    # use offsets to keep float conditioning: variables are O(1), kappa ~ 1.7e-3
    res = linprog([0, 0, 0, -1], A_ub=A, b_ub=bb, bounds=[(None, None)]*4, method='highs')
    assert res.status == 0, res.message
    v = res.x
    slack = [float(sum(F(x)*0 for x in co)) for _, co, _ in K]
    resid = [sum(float(co[i])*v[i] for i in range(4)) + float(k) for _, co, k in K]
    order = sorted(range(len(K)), key=lambda i: abs(resid[i]))
    # pick 4 tight, linearly independent constraints; solve exactly; verify feasibility exactly
    best = None
    for comb in itertools.combinations([i for i in order if abs(resid[i])<1e-9], 4):
        M = [list(K[i][1]) for i in comb]; rhs = [-K[i][2] for i in comb]
        # exact Gaussian elimination
        n = 4; Mx = [row[:] + [rhs[j]] for j, row in enumerate(M)]
        ok = True
        for col in range(n):
            piv = next((r for r in range(col, n) if Mx[r][col] != 0), None)
            if piv is None: ok = False; break
            Mx[col], Mx[piv] = Mx[piv], Mx[col]
            for r in range(n):
                if r != col and Mx[r][col] != 0:
                    f = Mx[r][col] / Mx[col][col]
                    Mx[r] = [x - f*y for x, y in zip(Mx[r], Mx[col])]
        if not ok: continue
        sol = [Mx[i][n] / Mx[i][i] for i in range(n)]
        if all(sum(co[i]*sol[i] for i in range(4)) + k >= 0 for _, co, k in K):
            if best is None or sol[3] > best[0][3]: best = (sol, [K[i][0] for i in comb])
    if best is not None:
        assert abs(float(best[0][3]) - v[3]) < 1e-12, (float(best[0][3]), v[3])
        return best
    raise SystemExit(("no exact vertex", [(K[i][0], resid[i]) for i in order[:10]], list(v)))
base = dict(c=F(3, 32), halo=48, core=128, union='radix', guard='fine')
cfgs = [
 ('FINAL (PR2044 rows, slacks->0)', base),
 ('- fine primes (gauge-free guard with theta~1/d), Strassen c=3/32', dict(base, guard='theta_d')),
 ('- fine primes, structured boundary c=1/2 (=PR1758 config, free eps)', dict(base, guard='theta_d', c=F(1, 2))),
 ('- fine primes - gauge-free (linear guard), Strassen', dict(base, guard='linear')),
 ('- radix union (old row 1-eps-beta-1/2+xi)', dict(base, union='old')),
 ('- per-visit tag (old 1-2beta)', dict(base, tag2=True)),
 ('+ old per-axis scan row 1-eps', dict(base, old_scan=True)),
 ('PR1673 strassen rows (linear guard, old union + full key, halo 3/core 8)', dict(c=F(3, 32), halo=3, core=8, union='old', full_key=True, guard='linear')),
 ('PR1673 + radix union (pkg4-like, Strassen)', dict(c=F(3, 32), halo=3, core=8, union='radix', guard='linear', tag2=True)),
 ('PR1673 + radix union + gauge-free (theta~1/d), Strassen', dict(c=F(3, 32), halo=3, core=8, union='radix', guard='theta_d')),
 ('PR1673 + radix union + gauge-free + structured c=1/2', dict(c=F(1, 2), halo=3, core=8, union='radix', guard='theta_d')),
]
out = []
fin = None
for name, cfg in cfgs:
    sol, act = solve(cfg)
    k = sol[3]
    if fin is None: fin = k
    rec = dict(config=name, sup_kappa='%.25f' % float(k) if False else str(+__import__('decimal').Decimal(k.numerator) / __import__('decimal').Decimal(k.denominator))[:32],
               q_minus_sup='%.10e' % float(q - k), final_minus_sup='%.10e' % float(fin - k), sup_minus_supplier='%.10e' % float(k - online),
               eps='%.12f' % float(sol[0]), beta='%.12f' % float(sol[1]), xi='%.12f' % float(sol[2]), active=act)
    out.append(rec); print(json.dumps(rec)); sys.stdout.flush()
json.dump(out, open(sys.argv[2], 'w'), indent=1)
