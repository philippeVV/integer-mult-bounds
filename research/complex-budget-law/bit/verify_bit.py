#!/usr/bin/env python3
"""Bit-side budget law for the selected PR589 (on PR555) bit supplier: a* ~ nu/Lambda.

usage: python3 -I -B verify_bit.py [--root REPO_ROOT] [--fast]
REPO_ROOT holds origin/main 0fdd7d82 files (default: ./main next to this script, a `git archive` export).
Standard library only. Own outward-rounded rational interval engine (2^-256), independent of the pinned
moment.py; the pinned engine is used only as a cross-check of the record endpoint.
"""
import gzip, hashlib, importlib.util, json, sys
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from math import isqrt, log, exp as fexp
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = Path(sys.argv[sys.argv.index('--root') + 1]) if '--root' in sys.argv else HERE.parents[2]
FAST = '--fast' in sys.argv
PINS = {
    'research/e8-bit-p111-tuned/tune/bitE8mut-p111-tuned.json.gz': 'e3e301e62bc2b41dbf1ceac0acf3d35824ebeab1848c5a2142d4577b109da52f',
    'research/e8-bit-p111-tuned/tune/EXPECTED.json': 'b9f0e5ca038abe837b0007c8fe006f366a6f02009244ac57705a8a200e2301f4',
    'research/e8-bit-p111-tuned/tune/wlib.py': 'a37fdf3cdf5f62945c2a90193401c0e9edc9dd3cda8d8cb36a3b9c329ebe4dbc',
    'research/e8-bit-p111-tuned/tune/fastc.py': 'ae51f1fbb22b8c8dac8aab787c798361f6aa5d4e1fefb29e43d0197899d4eab5',
    'research/e8-bit-carrier-skeleton/code/hub/e8-r740-replica32/code/engines/moment.py': '81cbc864e106418fc7854322746f72180cb8110ff39f74450a3b8d55a82c1482',
    'research/community-frontier1094-audit/ASSEMBLY.json': '1122cb3a575f8c788345fc1bf4c108ac7c3250fb81a5988671f5eb04b3ab3045',
}
for rel, h in PINS.items():
    assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == h, rel


def load(rel, name):
    sp = importlib.util.spec_from_file_location(name, ROOT / rel); mod = importlib.util.module_from_spec(sp); sp.loader.exec_module(mod); return mod


out = {}
f = lambda x, n=12: f'{float(x):.{n}g}'

# ------------------------------------------------------------------ interval engine (own)
PREC = 256
down = lambda x: Q((x.numerator << PREC) // x.denominator, 1 << PREC)
up = lambda x: -down(-x)


@lru_cache(None)
def atanh(t):
    s, p, k = Q(0), t, 0
    while True:
        s += p / (2 * k + 1); p *= t * t; k += 1
        tail = p / ((2 * k + 1) * (1 - t * t))
        if tail < Q(1, 1 << (PREC + 8)):
            return down(s), up(s + tail)


@lru_cache(None)
def ln(k):  # integer k >= 1
    if k == 1:
        return Q(0), Q(0)
    lo, hi = ln(k - 1); a, b = atanh(Q(1, 2 * k - 1))
    return down(lo + 2 * a), up(hi + 2 * b)


def expi(lo, hi):
    def series(y, upper):
        s, t, k = Q(1), Q(1), 0
        while True:
            k += 1; t = t * y / k; s += t
            if t * y < Q(1, 1 << (PREC + 8)):
                return up(s + 2 * t * y / (k + 1)) if upper else down(s)
    return series(lo, False), series(hi, True)


@lru_cache(None)
def Lr(r, m):  # enclosure of ln(m/r)
    lm, um = ln(m); lr, ur = ln(r)
    return lm - ur, um - lr


def psi(P, m, W, a, surcharge=True):
    """Enclosure of the bit moment (1/(mW)) [sum r n_r e^{a L_r} + eps m^a], eps = 32 m^2 N 1e-16 (pinned moment.py, bad=True)."""
    lo = hi = Q(0)
    for r, n in P.items():
        l, u = Lr(r, m); e = expi(a * l, a * u)
        lo += r * Q(n) * e[0]; hi += r * Q(n) * e[1]
    if surcharge:
        N = sum(Q(n) for n in P.values()); eps = Q(32 * m * m) * N / 10 ** 16
        l, u = ln(m); e = expi(a * l, a * u); lo += eps * e[0]; hi += eps * e[1]
    return down(lo / (m * W)), up(hi / (m * W))


def froot(P, m, W, surcharge=True):
    N = sum(float(n) for n in P.values()); eps = 32 * m * m * N / 1e16 if surcharge else 0
    g = lambda a: sum(r * float(n) * (m / r) ** a for r, n in P.items()) + eps * m ** a - m * W
    lo, hi = 0.0, 0.2
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if g(mid) < 0 else (lo, mid)
    return lo


GRID = Q(1, 10 ** (9 if FAST else 12))


def root(P, m=64, W=None, surcharge=True, grid=GRID):
    """Certified grid bracket [lo, lo+grid]: psi(lo) < 1 < psi(lo+grid)."""
    W = W or W0
    a = Q(int(froot(P, m, W, surcharge) / float(grid))) * grid
    for _ in range(50):
        if psi(P, m, W, a, surcharge)[1] < 1:
            break
        a -= grid
    while psi(P, m, W, a + grid, surcharge)[0] <= 1:
        a += grid
    assert psi(P, m, W, a, surcharge)[1] < 1 < psi(P, m, W, a + grid, surcharge)[0]
    return a


def passes(P, a, m=64, W=None):
    return psi(P, m, W or W0, a)[1] < 1


def logmass(P, m=64):
    return sum(r * Q(n) * Lr(r, m)[0] for r, n in P.items())


def plus(*hs, scale=None):
    c = defaultdict(Q)
    for h in hs:
        for r, n in h.items():
            c[r] += Q(n)
    return {r: n for r, n in sorted(c.items()) if n}


def times(h, k):
    return {r: Q(n) * k for r, n in h.items()}


# ------------------------------------------------------------------ 1. the word and the coverage law
wlib = load('research/e8-bit-p111-tuned/tune/wlib.py', 'wlib')
w = wlib.load(ROOT / 'research/e8-bit-p111-tuned/tune/bitE8mut-p111-tuned.json.gz')
EXP = json.loads((ROOT / 'research/e8-bit-p111-tuned/tune/EXPECTED.json').read_text())
ASM = json.loads((ROOT / 'research/community-frontier1094-audit/ASSEMBLY.json').read_text())
F = wlib.Frames(w); dim = F.dim
V, R, H = w['v'], w['R'], w['h']
m, L, K, T = 64, 9, 4, 32                      # m = 1 + 7L, K = (L-1)/2, 32 replicas
NI = L * T                                     # 288 local invocations per call
W0 = 2 * K * V * T                             # live families (all helpers hosted, no banks)
CST = 9 * 7
assert m == 1 + 7 * L and W0 == 28416 == EXP['nine']['bit']['live']

# own replay: per-register climb lists (same climb rule as wlib.replay), copies = retained totals' dims at scatter
N = 2 * V + R
cur = list(w['start']); climbs = defaultdict(list)


def climb(r, fr):
    a = cur[r]
    if a == fr:
        return
    assert F.le(a, fr)
    if dim[fr] > dim[a]:
        climbs[r].append(dim[fr] - dim[a])
    cur[r] = fr


for sec in ('D', 'A'):
    for fr, t, ss in w[sec]:
        for r in [t, *ss]:
            climb(r, fr)
copies = [dim[cur[q]] for q in w['ret']]
for sec in ('B', 'C'):
    for fr, t, ss in w[sec]:
        for r in [t, *ss]:
            climb(r, fr)
for r, fr in enumerate(w['final']):
    climb(r, fr)
cls = lambda r: 'source' if r < V else 'target' if r < 2 * V else ('total' if r in w['ret'] else 'helper')
HC = {c: Counter() for c in ('source', 'target', 'helper', 'total')}
ncl = {c: Counter() for c in HC}
for r in range(N):
    for d in climbs[r]:
        HC[cls(r)][d] += 1
    ncl[cls(r)][len(climbs[r])] += 1
HC['copy'] = Counter(copies)
local = plus(*HC.values())
assert local == {int(k): v for k, v in EXP['word']['local_histogram'].items()} == wlib.paid(w, F)
assert copies == [7] * 9 and sum(copies) == CST == EXP['word']['decoder_constant']
assert all(sum(climbs[r]) == 7 for r in range(2 * V)) and all(sum(climbs[r]) == 8 for r in range(2 * V, N))
lmass = sum(r * n for r, n in local.items())
assert lmass == 8 * R + 2 * V * 7 + CST == 5129

# slot model of the nine-stage hosting (PR1023 / bit589-width-packing): per pair, idle windows of 7*Delta,
# Delta in {1,1,1,2,2,4,4,6,6}; width-8 residual fits floor(7 Delta/8) = Delta-1 per window.
WIN = {7: 3, 14: 2, 28: 2, 42: 2}                         # windows per family, per call
slots = {c: 2 * V * k for c, k in WIN.items()}           # per replica: {7:666,14:444,28:444,42:444}
cap8 = {c: c // 8 for c in WIN}
per_rep = L * R                                          # 3951 width-8 residuals per replica
capacity = sum(slots[c] * cap8[c] for c in WIN)
assert capacity == 4 * V * L and per_rep <= capacity     # hosting wall R <= 4v
slack = capacity - per_rep                               # 45 residual slots unused
assert slack % cap8[42] == 0
rem_rep = Counter({7: slots[7], 14 - 8: slots[14], 28 - 24: slots[28], 42 - 40: slots[42] - slack // 5, 42: slack // 5})
REM = {r: n * T for r, n in rem_rep.items()}
fastc = (ROOT / 'research/e8-bit-p111-tuned/tune/fastc.py').read_text()
assert 'REM={2:13920,4:14208,6:14208,7:21312,42:288}' in fastc and REM == {2: 13920, 4: 14208, 6: 14208, 7: 21312, 42: 288}
TEMPLATE_IDLE = {c: k * 2 * V * T for c, k in WIN.items()}  # helper-free idle windows (R = 0)
assert sum(r * n for r, n in TEMPLATE_IDLE.items()) == NI * 8 * R + sum(r * n for r, n in REM.items())

G = plus(times(local, NI), REM)
assert G == {int(k): v for k, v in EXP['nine']['profile'].items()}
mass = sum(r * n for r, n in G.items())
D = m * W0 - mass
assert D == T * (2 * K * V - L * CST) == W0 - NI * CST == 10272
nu = Q(D, W0)
assert nu == Q(321, 888) == 1 - Q(L, 2 * K) * Q(CST, V)
out['coverage'] = dict(m=m, W=W0, invocations_per_call=NI, mass=mass, deficit=D, nu=str(nu), gamma=f'{CST}/{V}',
                       law='nu = u(1 - (L/(2K)) gamma) = 1 - (9/8)(63/111)', per_family_mass='63 + 63*288/W',
                       slot_model=dict(windows_per_family=WIN, width8_capacity_per_replica=capacity, residuals_per_replica=per_rep,
                                       unused_slots=slack, REM_per_replica=dict(sorted(rem_rep.items()))))

# ------------------------------------------------------------------ 2. tilt identity, exact root, sandwich
c_rec = Q(ASM['bit_lower'])
assert c_rec == Q(EXP['nine']['bit_lower']) == Q(1733427039589844234749, 10 ** 24)
up_rec = c_rec + Q(1, 10 ** 24)
lo_v, hi_v = psi(G, m, W0, c_rec), psi(G, m, W0, up_rec)
assert lo_v[1] < 1 < hi_v[0]                               # own engine: accept c, reject c + 1e-24
Mo = load('research/e8-bit-carrier-skeleton/code/hub/e8-r740-replica32/code/engines/moment.py', 'mom')
cert = Mo.certify({r: int(n) for r, n in G.items()}, m, W0, True)
assert Q(cert['lower']) == c_rec                           # pinned engine proposes the same endpoint
lam_lo = logmass(G); lam_hi = sum(r * n * Lr(r, m)[1] for r, n in G.items())
Ncalls = sum(G.values()); eps = Q(32 * m * m * Ncalls, 10 ** 16)
ellbar = lam_lo / mass
Lmax = Lr(1, m)
delta = Q(D) / mass
U_lin = Q(D) / lam_lo                                       # a* <= D/(W Lambda) = nu/Lambda
U_jen = log(1 + float(delta)) / float(ellbar)     # Jensen (float; checked below to bracket)
L_chord = log(1 + float(delta) * float(Lmax[0]) / float(ellbar)) / float(Lmax[1])
L_exp = float(U_lin) * fexp(-float(U_lin) * float(Lmax[1]))
lam2 = sum(r * n * float(Lr(r, m)[0]) ** 2 for r, n in G.items())
second = float(U_lin) * (1 - float(U_lin) * lam2 / (2 * float(lam_lo)))
cf = float(c_rec)
assert L_exp < L_chord < cf < U_jen < float(U_lin)
# exact form of the identity at the endpoint bracket: sum r n (e^{aL}-1) + eps m^a  crosses D
out['theorem1'] = dict(
    identity='Psi(a)=1  <=>  sum_r r n_r (e^{a ln(m/r)} - 1) + eps m^a = D,  eps = 32 m^2 N/1e16 = %s' % f(eps, 6),
    record_endpoint=str(c_rec), own_engine_accepts_c_rejects_next=True, pinned_engine_same_endpoint=True,
    Lambda_per_family=f(lam_lo / W0), ellbar=f(ellbar), r_geo=f(m * fexp(-float(ellbar))), delta=f(delta),
    upper_linear_nu_over_Lambda=f(U_lin), upper_jensen=f(U_jen), lower_chord=f(L_chord), lower_exp=f(L_exp),
    second_order=f(second), rel_gap_linear=f((float(U_lin) - cf) / cf, 4), rel_gap_second=f((second - cf) / cf, 4),
    surcharge_shift_of_root=f(float(eps) * 64 ** cf / float(lam_lo), 3))

# complex side, same law, for side-by-side
CPLX = {1: 349604, 2: 60192, 3: 62585, 4: 20448, 5: 11808, 7: 7, 8: 1152, 48: 4866}
cl = sum(r * n * Lr(r, 73)[0] for r, n in CPLX.items())
out['complex_side_same_law'] = dict(nu='5/14', Lambda=f(cl / 14336), nu_over_Lambda=f(Q(5, 14) / (cl / 14336)), a='0.0016754627')

# eleven-stage control: m = 78, W = 8880 (8 replicas), D = T(2Kv - L cst)
P11 = {int(k): v for k, v in EXP['eleven']['profile'].items()}
m11, W11 = 78, EXP['eleven']['bit']['live']
D11 = m11 * W11 - sum(r * n for r, n in P11.items())
assert D11 == 8 * (2 * 5 * V - 11 * CST)
c11 = Q(EXP['eleven']['bit_lower'])
assert psi(P11, m11, W11, c11)[1] < 1 < psi(P11, m11, W11, c11 + Q(1, 10 ** 24))[0]
out['eleven_stage_control'] = dict(nu=str(Q(D11, W11)), Lambda=f(logmass(P11, m11) / W11), a=str(c11),
                                   note='more stages: nu up 4% (11/10 vs 9/8) but Lambda up 6.7% (m=78, idle per family 7(L-2.2))')

# ------------------------------------------------------------------ 3. decomposition of Lambda (per family)
NIQ = Q(NI)
parts = {}
idle_parts = {'idle 7-windows (unhosted)': {7: REM[7]}, 'remainders in hosted windows (6/4/2)': {6: REM[6], 4: REM[4], 2: REM[2]},
              'empty 42-windows': {42: REM[42]}}
groups = [('data: sources (111 x 7)', times(HC['source'], NIQ)), ('data: targets (111 x 7)', times(HC['target'], NIQ)),
          ('helpers (430 x 8)', times(HC['helper'], NIQ)), ('retained totals (9 x 8)', times(HC['total'], NIQ)),
          ('decoder copies (9 x rank 7)', times(HC['copy'], NIQ))] + list(idle_parts.items())
tot = Q(0)
for name, h in groups:
    lm_ = logmass(h); tot += lm_
    parts[name] = dict(mass_per_family=f(sum(r * n for r, n in h.items()) / W0, 5), log_per_family=f(lm_ / W0, 5),
                       share=f(lm_ / lam_lo, 4), hist=dict(sorted((r, f(n / NI, 6)) for r, n in h.items())))
assert tot == lam_lo
tmpl = logmass(TEMPLATE_IDLE)
hosted = logmass(times(plus(HC['helper'], HC['total']), NIQ)) + logmass(REM)
out['decomposition'] = parts
out['hosting_fragmentation'] = dict(
    helper_free_template_idle_log_per_family=f(tmpl / W0, 5), hosted_helpers_plus_remainder_log_per_family=f(hosted / W0, 5),
    fragmentation_cost_share_of_Lambda=f((hosted - tmpl) / lam_lo, 4),
    note='helpers are idle-displacing: their 8R climbs are carved out of the template idle windows, so they add no mass and no '
         'deficit, only log-price (fragmentation of rank-7..42 windows into rank-1/2 climbs)')
lam_data = logmass(times(plus(HC['source'], HC['target']), NIQ)); lam_copy = logmass(times(dict(HC['copy']), NIQ))
assert lam_data + lam_copy + hosted == lam_lo
phi = (hosted - tmpl) / lam_lo
out['Lambda_partition_per_family'] = dict(data=f(lam_data / W0, 5), copies=f(lam_copy / W0, 5), template_windows=f(tmpl / W0, 5),
                                          helper_fragmentation=f((hosted - tmpl) / W0, 5), total=f(lam_lo / W0, 6),
                                          R0_first_order_multiplier_1_over_1_minus_phi=f(1 / (1 - phi), 5))
out['exchange_rates_first_order'] = dict(
    copy_unit_per_invocation_pct_of_a=f(100 * Q(NI) / (W0 * nu), 5),
    plane_port_per_invocation_pct_of_a=f(100 * Q(8, 9) * Q(NI) / (W0 * nu), 5),
    helper_index_pct_of_a=f(100 * (hosted - tmpl) / R / lam_lo, 4),
    copy_unit_in_helper_indices=f((Q(NI) / (W0 * nu)) / ((hosted - tmpl) / R / lam_lo), 4))
out['climb_counts'] = {c: dict(sorted(v.items())) for c, v in ncl.items()}

# ------------------------------------------------------------------ 4. lever ceilings (exact grid brackets)
data = plus(HC['source'], HC['target']); helpers = plus(HC['helper'], HC['total']); cop = dict(HC['copy'])


def prof(dat=data, hel=helpers, cp=cop, idle=REM):
    return plus(times(dat, NIQ), times(hel, NIQ), times(cp, NIQ), idle)


def extreme(cls_list):
    h = Counter()
    for r in range(N):
        if cls(r) in cls_list:
            c, wd = len(climbs[r]), sum(climbs[r])
            h[1] += c - 1; h[wd - c + 1] += 1
    return {k: v for k, v in h.items() if v}


def idle_to(rank, extra_mass):
    return {rank: Q(extra_mass) / rank}


hmass = sum(r * n for r, n in helpers.items()) * NI
base = root(G)
rows = [
    ('baseline (selected profile, slot-exact hosting, surcharge)', G),
    ('idle packing: all REM at rank 45 (free concentration; PR1023 proves width-8 slot hosting is already optimal)',
     plus(times(local, NIQ), idle_to(45, sum(r * n for r, n in REM.items())))),
    ('re-ranking with fixed climb counts, data+helpers ((c-1) unit + one big; Lemma 3 extreme)',
     prof(dat=extreme({'source', 'target'}), hel=extreme({'helper', 'total'}))),
    ('one-shot data: every source/target one rank-7 child per invocation', prof(dat={7: 2 * V})),
    ('one-shot helpers: every helper one rank-8 child', prof(hel={8: R})),
    ('one-shot data + helpers', prof(dat={7: 2 * V}, hel={8: R})),
    ('all 439 helpers width 7 (slot optimum {42:340,21:1}/replica; rank-1 child removed) [PR1023 control +23.41%]',
     prof(hel=plus(helpers, {1: -R}), idle={42: 340 * T, 21: T})),
    ('helpers removed (R=0): template windows {7,14,28,42} idle', prof(hel={}, idle=TEMPLATE_IDLE)),
    ('helpers removed + data one-shot', prof(dat={7: 2 * V}, hel={}, idle=TEMPLATE_IDLE)),
    ('decoder free (cst=0, nu=1)', prof(cp={})),
    ('decoder at the bit label floor gamma=9/17 (cst/v; nu=55/136), carried as rank-7 copies', None),
]
lever = []
for name, P in rows:
    if P is None:     # gamma = 9/17 at v = 111: cst = 999/17 carried as rank-7 copies (count cst/7)
        cst = Q(9, 17) * V
        P = prof(cp={7: cst / 7})
    a = base if P is G else root(P)
    lever.append(dict(lever=name, a=f(a), x=f(a / base, 4), kappa_cap=f(a / (1 + a))))
out['levers'] = lever
assert abs(float(lever[6]['a']) - 0.0021391622) < 2e-10 + 2 * float(GRID)  # bit589-width-packing width-7 row

# universal cap of the template: all paid mass at the halving cap 45 (row bound forces HD=2)
for cst in (CST, 0):
    Pu = {45: Q(63 * W0 + NI * cst, 45)}
    out[f'universal_cap_cst{cst}'] = f(root(Pu))

# ------------------------------------------------------------------ 5. requirements for target savings
TARGETS = {'a=0.0025': Q(25, 10000), 'a=0.0034': Q(34, 10000), 'complex kappa x2 (a=0.0033565)': None}
KAP = Q(836330098538459, 500000000000000000)
TARGETS['complex kappa x2 (a=0.0033565)'] = 2 * KAP / (1 - 2 * KAP)


def bisect(pred, lo, hi, it=40):
    if not pred(lo):
        return None
    for _ in range(it):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if pred(mid) else (lo, mid)
    return lo


def plane(cst2, shrink, ports=None):
    """Plane ports (ctoth #1587 / plane-port-copy-budgets) on `ports` of the v live ports (default all, u = 2):
    each upgraded port's data incidences climb 6 per stage; copies cst2 at rank 7; deficit (1 + ports/v) W - NI cst2.
    Idle = remaining mass: surplus at rank 42 (freed windows), shortfall by scaling REM proportionally."""
    k = Q(V if ports is None else ports)
    dat = plus(data, {1: -2 * k}) if shrink else data
    paid = sum(r * n for r, n in plus(times(dat, NIQ), times(helpers, NIQ)).items()) + NI * cst2
    idle_mass = (m - 1 - k / V) * W0 + NI * cst2 - paid
    rem_mass = sum(r * n for r, n in REM.items())
    idle = plus(REM, {42: (idle_mass - rem_mass) / 42}) if idle_mass >= rem_mass else times(REM, idle_mass / rem_mass)
    return plus(times(dat, NIQ), times(helpers, NIQ), {7: NIQ * cst2 / 7} if cst2 else {}, idle)


assert plane(Q(CST), False, 0) == G                       # zero upgraded ports reproduces the record profile


need = {}
for tname, A in TARGETS.items():
    rec = {'a_target': f(A), 'x_over_record': f(A / c_rec, 4), 'Lambda_needed_at_nu_321/888 (first order)': f(nu / A, 5)}
    cmax = bisect(lambda c: passes(prof(cp={7: Q(c) / 7} if c else {}), A), Q(0), Q(63))
    rec['max_cst_u1 (shipped word)'] = None if cmax is None else f(cmax, 5)
    cmax7 = bisect(lambda c: passes(prof(hel=plus(helpers, {1: -R}), cp={7: Q(c) / 7} if c else {}, idle={42: 340 * T, 21: T}), A), Q(0), Q(63))
    rec['max_cst_u1 (all helpers width 7)'] = None if cmax7 is None else f(cmax7, 5)
    t = bisect(lambda t: passes(prof(hel=times(helpers, t), idle=plus(REM, idle_to(42, (1 - t) * hmass))), A), Q(0), Q(1))
    rec['max_helper_climb_fraction (freed mass to rank 42)'] = None if t is None else f(t, 4)
    t0 = bisect(lambda t: passes(prof(hel=times(helpers, t), cp={}, idle=plus(REM, idle_to(42, (1 - t) * hmass))), A), Q(0), Q(1))
    rec['max_helper_climb_fraction at cst=0'] = None if t0 is None else f(t0, 4)
    lm_ = sum(r * n for r, n in plus(data, helpers).items())
    rec['min_uniform_local_rank (data+helpers, cst 63)'] = next((r for r in range(1, 46) if passes(prof(dat={r: Q(lm_, r)}, hel={}), A)), None)
    for sh in (False, True):
        c2 = bisect(lambda c: passes(plane(Q(c), sh), A), Q(0), Q(400))
        rec[f'max_cst2_u2 ({"shrink" if sh else "unchanged"})'] = None if c2 is None else f(c2, 5)
    for sh in (False, True):       # partial plane upgrade at unchanged copy charge: minimal number of upgraded ports
        k = bisect(lambda k: not passes(plane(Q(CST), sh, k), A), Q(0), Q(V))
        rec[f'min_upgraded_ports_at_dC0 ({"shrink" if sh else "unchanged"})'] = None if k is None else f(k, 5)
    need[tname] = rec
out['requirements'] = need
# combinations inside the line-port law (u = 1)
FLOOR = Q(9, 17) * V                                      # bit label-decoder ratio floor cst/v >= 9/17 (bit-orthogonal-pair-decoders)
w7 = dict(hel=plus(helpers, {1: -R}), idle={42: 340 * T, 21: T})
combos = [('width-7 helpers + decoder at the 9/17 floor', prof(cp={7: FLOOR / 7}, **w7)),
          ('one-shot data + one-shot helpers + 9/17 decoder', prof(dat={7: 2 * V}, hel={8: R}, cp={7: FLOOR / 7})),
          ('width-7 helpers + one-shot data + 9/17 decoder', prof(dat={7: 2 * V}, cp={7: FLOOR / 7}, **w7)),
          ('everything line-port: one-shot data, width-7 one-shot helpers, 9/17 decoder',
           prof(dat={7: 2 * V}, hel={7: R}, cp={7: FLOOR / 7}, idle={42: 340 * T, 21: T}))]
out['line_port_combinations'] = [dict(combo=n, a=f(root(P)), x=f(root(P) / base, 4)) for n, P in combos]

out['plane_port_ceilings'] = {f'cst2={c} {"shrink" if s else "unchanged"}': f(root(plane(Q(c), s))) for c in (0, 63) for s in (False, True)}

print(json.dumps(out, indent=1, default=str))
print('PASS_BIT_BUDGET_LAW')
