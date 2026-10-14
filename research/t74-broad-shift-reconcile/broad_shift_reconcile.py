#!/usr/bin/env python3
"""Independent exact reconciliation of the broad-shift correction with the inherited outer
(PR2044 head 322949a5, research/capped-fusion-t74-k554/outer/GAUGEFREE-INTERIORS.md sec. 2, 5;
inherited definitions from vendor/pr1673/.../CAPPED-GAUGE-AUDIT.md sec. 3 and FORWARD-SMALL-A-AUDIT.md sec. 2).
No producer imports. Part A: sympy identities (exact rational functions). Part B: exhaustive exact
Fraction checks over every phase cell / block / anchor on a finite parameter grid. b0 normalised to 1
(every exponent is b0 * polynomial), alpha^-1 constant omitted.
Usage: python broad_shift_reconcile.py [qmax] [cden]"""
import sys, json, hashlib, itertools
from fractions import Fraction as F
from pathlib import Path
import sympy as sp

def req(ok, why):
    if not ok: raise SystemExit('FAIL: ' + why)

R = {}
# ---------------- Part A: symbolic ----------------
x, y, b, th, c, m, s, w = sp.symbols('x y b theta c m s w', real=True)
rho = 1 + th; lam = 1 / rho
dl = sp.Symbol('delta', real=True)
EK = -b * (x - lam * y + dl) ** 2                       # broad kernel exponent (FORWARD-SMALL-A sec.2)
EG = lambda sh: -b * lam * (x - y - sh) ** 2              # shifted Toeplitz kernel exponent
z = lambda e: sp.simplify(sp.together(sp.expand(e))) == 0
# A1 Lemma 2.1 for every s
lem = -b * ((1 - lam) * x**2 - lam * (1 - lam) * y**2 + dl**2 - lam * s**2 + 2 * (dl + lam * s) * x - 2 * lam * (dl + s) * y)
req(z(EK - EG(s) - lem), 'Lemma 2.1 general-s identity')
# A2 inherited delta=(1-lambda)m-c and beta'_m=theta m - rho c  ==>  delta = beta'_m / rho
delta_def = (1 - lam) * m - c; betap = th * m - rho * c
req(z(delta_def - betap / rho), 'delta = beta_prime_m / rho')
R['A2'] = 'delta = beta\'_m/rho exactly; beta\'_m in [-1/2,1/2) gives |delta| <= 1/(2 rho) < 1/2 (inherited |delta|<=1/2 holds, strictly)'
# A3 unshrunk shift s1=-delta/lambda = -beta'_m
s1 = -dl / lam
sub1 = sp.expand(EK - EG(s1))
correct1 = -b * ((1 - lam) * x**2 - lam * (1 - lam) * y**2 + dl**2 - dl**2 / lam + 2 * dl * (1 - lam) * y)
req(z(sub1 - correct1), 'correct s=-delta/lambda bracket')
inherited_or_true = -b * ((1 - lam) * x**2 - lam * (1 - lam) * y**2 + dl**2 + dl**2 / lam - 2 * dl * (1 - lam) * y)
resid = sp.factor(sp.simplify(sub1 - inherited_or_true))
req(resid != 0, 'inherited or-True expression is wrong')
R['A3_inherited_or_True_residual'] = str(resid)
req(z(sp.diff(sub1, x).subs({x: 0, y: 0})), 'x-linear vanishes at s1')
R['A3'] = '|s1| = |beta\'_m| <= 1/2, equality at beta\'_m=-1/2 (attained, cell is half-open) -> shrink needed'
# symbol vanishes at |s|=1/2, w=pi: pairing u <-> 1-u in sum_u exp(-g(u-1/2)^2)(-1)^u
u = sp.Symbol('u', integer=True)
req(sp.simplify(((u - sp.Rational(1, 2))**2) - (((1 - u) - sp.Rational(1, 2))**2)) == 0 and sp.simplify((-1)**u + (-1)**(1 - u)) == 0,
    'phi(pi)=0 at s=1/2 by pairing')
# A4 shrunk shift used by the outer: s2=-delta(1-theta)
s2 = -dl * (1 - th)
sub2 = sp.expand(EK - EG(s2))
cx = sp.simplify(sp.diff(sub2, x).subs({x: 0, y: 0})); cy = sp.simplify(sp.diff(sub2, y).subs({x: 0, y: 0}))
req(z(cx - (-b) * 4 * dl * th / (1 + th)), 'x-linear 4 b delta theta/(1+theta)')
req(z(cy - (-b) * (-2 * dl * th / (1 + th))), 'y-linear -2 b delta theta/(1+theta)')
# with |delta|<=1/(2rho): |s2| <= (1-th)/(2(1+th)), 1-2|s2| >= 2th/(1+th) >= th  (0<th<=1)
gap = 1 - (1 - th) / (1 + th)
req(z(gap - 2 * th / (1 + th)), 'gap identity')
req(z((2 * th / (1 + th) - th) - th * (1 - th) / (1 + th)), 'gap >= theta for 0<theta<=1')
R['A4'] = ('s2=-delta(1-theta): x-lin 4b delta th/(1+th) (<= 2b th/(1+th)^2), y-lin 2b delta th/(1+th)^2 abs (<= b th/(1+th)^2); '
           '1-2|s2| >= 2th/(1+th) >= th. Outer uses s2, whose asserts in check_gaugefree.py stage A are live (not or-True).')
# A5 the diagonals: Dcol'(x) = -b[(1-lam)x^2 + 2(dl+lam s)x], Drow'(y) = -b[-lam(1-lam)y^2 - 2lam(dl+s)y + dl^2 - lam s^2]
logDcol = lambda sh: -b * ((1 - lam) * x**2 + 2 * (dl + lam * sh) * x)
logDrow = lambda sh: -b * (-lam * (1 - lam) * y**2 - 2 * lam * (dl + sh) * y + dl**2 - lam * sh**2)
req(z(EK - (logDrow(s2) + EG(s2) + logDcol(s2))), 'K = Drow G_s2 Dcol exact')
# narrow Lemma 2.2 (independent restatement)
a, be0, ti, h = sp.symbols('a beta0 t_i h', real=True)
bi = be0 + th * ti; bj = bi + th * h
EN = -a * ((rho * h + bi)**2 - bj**2)
req(z(EN - (-a * rho * h**2 - 2 * a * be0 * h + a * th * (ti**2 - (ti + h)**2))), 'Lemma 2.2')
req(z((-a * rho * h**2 - 2 * a * be0 * h) - (a * be0**2 / rho - a * rho * (h + be0 / rho)**2)), 'narrow shift s=-beta0/rho')
R['A_symbolic'] = 'PASS'

# ---------------- Part B: exhaustive exact grid ----------------
qmax = int(sys.argv[1]) if len(sys.argv) > 1 else 24
cden = int(sys.argv[2]) if len(sys.argv) > 2 else 12
half = F(1, 2)
cnt = dict(thetas=0, cells=0, blocks=0, entries=0)
max_abs_s1 = F(0); max_delta_ratio = F(0); min_gap_ratio = None
worst_total_lin = (F(0), None)   # max over blocks of (lin range Dcol + lin range Drow)/(theta*B)  vs sentence "<= 2 b0 theta B (total)"
worst_sec5 = (F(0), None)        # max over blocks/diagonals of range / (theta B^2/4 + 4 theta B)
thetas = sorted({F(p, q) for q in range(1, qmax + 1) for p in range(1, q + 1)})
for T in thetas:
    cnt['thetas'] += 1
    Rh = 1 + T; L = 1 / Rh; S1T = 1 - T
    for cn in range(-2 * cden, 2 * cden):
        C = F(cn, cden)
        # one phase cell: all integer k with beta'_k = T k - Rh C in [-1/2,1/2)
        import math
        k0 = math.ceil((Rh * C - half) / T); cell = []
        k = k0
        while T * k - Rh * C < half:
            if T * k - Rh * C >= -half: cell.append(k)
            k += 1
        if not cell: continue
        cnt['cells'] += 1
        n = len(cell)
        Bs = [1 << e for e in range(0, n.bit_length()) if (1 << e) <= n]
        for B in Bs:
            for st in range(0, n - B + 1):
                blk = cell[st:st + B]; M = blk[B // 2]          # anchor at centre
                cnt['blocks'] += 1
                D = (1 - L) * M - C
                bp = T * M - Rh * C
                req(D * Rh == bp, 'delta*rho == beta_prime')
                req(abs(D) <= half / Rh, '|delta| <= 1/(2rho)')
                max_delta_ratio = max(max_delta_ratio, abs(D) * 2 * Rh)
                S1 = -D / L; max_abs_s1 = max(max_abs_s1, abs(S1)); req(abs(S1) <= half, '|s1|<=1/2')
                S2 = -D * S1T
                g = 1 - 2 * abs(S2)
                req(g >= 2 * T / (1 + T) and g >= T and g > 0, 'gap 1-2|s2|')
                gr = g / T; min_gap_ratio = gr if min_gap_ratio is None else min(min_gap_ratio, gr)
                xs = [kk - M for kk in blk]
                dc = [-((1 - L) * X * X + 2 * (D + L * S2) * X) for X in xs]
                dr = [-(-L * (1 - L) * Y * Y - 2 * L * (D + S2) * Y + D * D - L * S2 * S2) for Y in xs]
                lc = [-(2 * (D + L * S2) * X) for X in xs]; lr = [(2 * L * (D + S2) * Y) for Y in xs]
                if B <= 16:   # entrywise exact factorisation on every entry of small blocks
                    for X in xs:
                        for Y in xs:
                            lhs = -(X - L * Y + D) ** 2
                            rhs = dr[xs.index(Y)] + (-L * (X - Y - S2) ** 2) + dc[xs.index(X)]
                            req(lhs == rhs, 'entrywise K = Drow G Dcol')
                            cnt['entries'] += 1
                tl = (max(lc) - min(lc)) + (max(lr) - min(lr))
                r = tl / (T * B)
                if r > worst_total_lin[0]: worst_total_lin = (r, dict(theta=str(T), c=str(C), B=B, m=M, delta=str(D), lin_range_total=str(tl)))
                bound5 = T * B * B / 4 + 4 * T * B
                for name, dd in (('Dcol', dc), ('Drow', dr)):
                    rr = (max(dd) - min(dd)) / bound5
                    req(rr <= 1, 'Section 5 per-diagonal range bound ' + name)
                    if rr > worst_sec5[0]: worst_sec5 = (rr, dict(diag=name, theta=str(T), c=str(C), B=B, m=M))
R['B_grid'] = dict(qmax=qmax, cden=cden, **cnt)
R['B_max_|s1|'] = str(max_abs_s1)
R['B_max_|delta|*2rho'] = str(max_delta_ratio)
R['B_min_(1-2|s2|)/theta'] = str(min_gap_ratio)
R['B_sec5_max_range_over_bound'] = [str(worst_sec5[0]), worst_sec5[1]]
R['B_lemma2.1_sentence_total_linear_range_over_theta_B'] = [str(worst_total_lin[0]), float(worst_total_lin[0]), worst_total_lin[1]]
R['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
R['sympy'] = sp.__version__
print(json.dumps(R, indent=1))
print('PASS_BROAD_SHIFT_RECONCILE')
