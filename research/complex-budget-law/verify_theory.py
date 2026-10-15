#!/usr/bin/env python3
"""Verify the deficit/log-mass law for the complex-side budget on the selected #1094 record.

usage: python3 -I -B verify_theory.py [--fast]
Inputs (sha256-pinned, from this repository): research/complex-local-word-ceiling/ceiling.py
and research/community-frontier1094-audit/selected-replay.tar.gz.
Standard library only; exact rational interval arithmetic via ceiling.py (outward rounding 2^-256).
"""
import gzip
import hashlib
import io
import json
import sys
import tarfile
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'complex-local-word-ceiling'))
import ceiling as C  # noqa: E402

PINS = {
    '../complex-local-word-ceiling/ceiling.py': '86a1432870b3b21d8454106ec9162fddc0471472d9214c590a1420f2426cedfc',
    '../community-frontier1094-audit/selected-replay.tar.gz': 'f1da80610bf33c58527f18ad63b4350d0b864e1d4486141ee17ac9bed4d529d6',
}
for rel, h in PINS.items():
    assert hashlib.sha256((HERE / rel).read_bytes()).hexdigest() == h, rel

M73, FAM, ST, CEN = 73, 896, 9, 64
W, NCALL = 14336, 144
GLOBAL = {1: 349604, 2: 60192, 3: 62585, 4: 20448, 5: 11808, 7: 7, 8: 1152, 48: 4866}
LOCAL = {1: 2412, 2: 418, 3: 432, 4: 142, 5: 82, 8: 8}
KAPPA = Q(836330098538459, 500000000000000000)
A_SEL_LO = Q(67018507072324140917, 40000000000000000000000)   # certified complex endpoint (accepted)
A_SEL_REJ = Q(837731338404051761463, 500000000000000000000000)  # next grid point (rejected)
FAST = '--fast' in sys.argv
GRID = Q(1, 10 ** (9 if FAST else 12))
out = {}


def f(x, n=12):
    return f'{float(x):.{n}g}'


def lnq(x):
    """Enclosure of ln x for rational x >= 1 via ln p - ln q for integers is too slow; use 2 atanh((x-1)/(x+1))."""
    x = Q(x)
    if x == 1:
        return Q(0), Q(0)
    t = (x - 1) / (x + 1)
    if t <= Q(1, 3):
        lo, hi = C.atanh(t)
        return C.down(2 * lo), C.up(2 * hi)
    raise ValueError('lnq range')


def L(r, m=M73):
    lm, um = C.ln(m)
    lr, ur = C.ln(r)
    return lm - ur, um - lr


def moment(hist, Wt, a, m=M73):
    lo = hi = Q(0)
    for r, n in hist.items():
        c = C.child_cost(r, m, a)
        lo, hi = lo + n * c[0], hi + n * c[1]
    return lo / Wt, hi / Wt


def root(hist, Wt, grid=Q(1, 10 ** 15), m=M73):
    lo, hi = Q(0), Q(1, 20)
    while hi - lo > grid:
        mid = Q(int((lo + hi) / 2 / grid)) * grid
        if mid <= lo:
            mid = lo + grid
        v = moment(hist, Wt, mid, m)
        if v[1] < 1:
            lo = mid
        elif v[0] > 1:
            hi = mid
        else:
            raise ValueError('enclosure too wide')
    return lo, hi


# ---------------------------------------------------------------- 1. coverage invariants
mass = sum(r * n for r, n in GLOBAL.items())
D = M73 * W - mass
assert mass == 72 * W + CEN * NCALL == 1041408 and D == 5120
nu = Q(D, W)
assert nu == Q(5, 14) == 1 - Q(9, 8) * Q(64, 112)          # nu = u(1 - (s/2P) * gamma)
assert Q(W, NCALL) == Q(FAM, ST)
out['coverage'] = dict(mass=mass, deficit=D, nu=str(nu), gamma='4/7', s_over_2P='9/8')

# ---------------------------------------------------------------- 2. theorem 1 (exact) on the selected profile
lam_lo = lam_hi = Q(0)
lam2 = 0.0
for r, n in GLOBAL.items():
    a, b = L(r)
    lam_lo += n * r * a
    lam_hi += n * r * b
    lam2 += n * r * float(a) ** 2
Lmax_lo, Lmax_hi = L(1)
U1 = (Q(D) / lam_lo, Q(D) / lam_hi)                     # a* <= D / Lambda
x = Q(D) / lam_lo
e = C.exp(Q(0), x * Lmax_hi)                              # a* >= (D/Lambda) m^{-D/Lambda}
L1 = (Q(D) / lam_hi) / e[1]
delta = Q(D, mass)
ld = lnq(1 + delta)
ellbar = (lam_lo / mass, lam_hi / mass)
U2 = ld[1] / ellbar[0]                                    # Jensen: a* <= ln(1+delta)/ellbar
# chord: a* >= ln(1 + delta Lmax/ellbar) / Lmax
y = delta * Lmax_lo / ellbar[1]
L2 = lnq(1 + y)[0] / Lmax_hi if y <= Q(1, 2) else None
r_lo, r_hi = root(GLOBAL, W)
assert r_lo <= A_SEL_REJ and r_hi >= A_SEL_LO
assert L1 <= r_lo and L2 <= r_lo and r_hi <= U2 <= U1[0]
lin = float(D) / float(lam_lo)
second = lin * (1 - lin * lam2 / (2 * float(lam_lo)))
out['theorem1'] = dict(
    Lambda_per_family=f(lam_lo / W), ellbar=f(ellbar[0]), r_geo=f(M73 / 2.718281828459045 ** float(ellbar[0])),
    delta=f(delta), root=[f(r_lo, 15), f(r_hi, 15)], selected_endpoint=f(A_SEL_LO, 15),
    upper_linear_D_over_Lambda=f(U1[0]), upper_jensen=f(U2), lower_chord=f(L2), lower_exp=f(L1),
    second_order_estimate=f(second), rel_gap_upper=f((U1[0] - r_lo) / r_lo, 4), rel_gap_second=f((second - float(r_lo)) / float(r_lo), 4))

# ---------------------------------------------------------------- 3. per-call log-mass decomposition from the shipped local word
tf = tarfile.open(HERE / '../community-frontier1094-audit/selected-replay.tar.gz')
member = 'candidate/strict/global/base/providers/gauge44/LOCAL-PROGRAMS.json.gz'
prog = json.load(gzip.open(io.BytesIO(tf.extractfile(member).read())))
dec = {}
widths = Counter()
for side in ('forward', 'backward'):
    H = {'data': Counter(), 'helper': Counter(), 'centre': Counter()}
    per_reg = Counter()
    nclimb = Counter()
    seg_touched = {}
    untouched = 0
    for ev in prog[side]:
        if ev[0] == 'climb':
            H['data' if ev[2] < 224 else 'helper'][ev[5]] += 1
            per_reg[ev[2]] += ev[5]
            nclimb[ev[2]] += 1
            if not seg_touched.get(ev[2], True):
                untouched += 1
            seg_touched[ev[2]] = False
        elif ev[0] == 'center':
            H['centre'][ev[4]] += 1
        elif ev[0] == 'scalar':
            for x in ev[1][1:]:
                if isinstance(x, list) and x[0] == 'L' and x[1] in seg_touched:
                    seg_touched[x[1]] = True
    assert untouched == 0  # every climb after a register's first is preceded by a scalar touch on it
    tot = H['data'] + H['helper'] + H['centre']
    assert dict(tot) == LOCAL, side
    assert all(per_reg[r] == 8 for r in range(224))
    dec[side] = H
    widths = Counter(v for r, v in per_reg.items() if r >= 224)
    widths_by_reg = dict(per_reg)
assert dec['forward'] == dec['backward']
H = dec['forward']
assert dict(widths) == {9: 240, 8: 122, 7: 30, 5: 32, 4: 56}
idle = dict((Counter(GLOBAL) - Counter({r: n * NCALL for r, n in LOCAL.items()})))
parts = {}
for k, h in list(H.items()) + [('idle(actual)', Counter({r: Q(n, NCALL) for r, n in idle.items()}))]:
    mlo = sum(n * r * L(r)[0] for r, n in h.items())
    parts[k] = dict(mass=f(sum(r * n for r, n in h.items())), logmass=f(mlo, 6), share=f(mlo * NCALL / lam_lo, 4),
                    hist={r: f(n, 6) for r, n in sorted(h.items())})
out['decomposition_per_call'] = parts
out['touch_separated_climbs'] = sum(nclimb.values())

# climb-count floor: keep each register's number of climbs c and width w; the cheapest ranks are
# (c-1) unit steps plus one step of w-c+1 (r ln(m/r) is concave), i.e. coarsening without new frames
floor_hist = Counter()
for r, c in nclimb.items():
    floor_hist[1] += c - 1
    floor_hist[per_reg[r] - c + 1] += 1
floor_hist = {k: v for k, v in floor_hist.items() if v}

# ---------------------------------------------------------------- 4. lever ceilings (ceiling.py relaxation, idle at R=51)
def cap(local, centre=CEN, unpaid=1, families=FAM):
    return C.ceiling(local, M73, families, ST, centre, grid=GRID, unpaid=unpaid)[1]


def plus(*hs):
    c = Counter()
    for h in hs:
        for r, n in h.items():
            c[r] += n
    return {r: n for r, n in c.items() if n}


data, helper, centre = dict(H['data']), dict(H['helper']), dict(H['centre'])
data1 = {8: 224}
helper1 = dict(widths)
base = cap(LOCAL)
rows = [
    ('baseline (selected word, ideal idle)', base),
    ('frame coarsening: every register climbs in ONE child (data+helpers)', cap(plus(data1, helper1, centre))),
    ('same climb counts, best ranks ((c-1) unit steps + one big)', cap(plus(floor_hist, centre))),
    ('  data registers only one-shot', cap(plus(data1, helper, centre))),
    ('  helper registers only one-shot', cap(plus(data, helper1, centre))),
    ('local-dimension cap: all non-centre local mass at rank h=9', cap(plus({9: Q(5522, 9)}, centre))),
    ('helper climbs removed (helper mass becomes idle)', cap(plus(data, centre))),
    ('helper climbs removed + data one-shot', cap(plus(data1, centre))),
    ('centre/decoder free, u=1 (C=0)', cap(plus(data, helper), centre=0)),
    ('plane ports u=2, decoder free (C=0)', cap(plus(data, helper), centre=0, unpaid=2)),
    ('plane ports u=2, C=128 as rank-2 terms (gamma=4/7)', cap(plus(data, helper, {2: 64}), centre=128, unpaid=2)),
    ('plane ports u=2, C=128 + helper one-shot', cap(plus(data, helper1, {2: 64}), centre=128, unpaid=2)),
]
out['levers'] = [dict(lever=n, a_cap=f(a), x=f(a / base, 4), kappa_cap=f(a / (1 + a))) for n, a in rows]

# ---------------------------------------------------------------- 5. what each lever needs for kappa x2 and kappa = 0.01
T2 = 2 * KAPPA / (1 - 2 * KAPPA)
T01 = Q(1, 99)


def ok(local, a, centre=CEN, unpaid=1):
    try:
        v = C.relax(local, M73, FAM, ST, centre, a, unpaid=unpaid)
    except ValueError:
        return False
    return v[1] < 1


def bisect(pred, lo, hi, it=40):
    """largest x in [lo,hi] with pred(x) True, assuming monotone True->False."""
    if not pred(lo):
        return None
    for _ in range(it):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if pred(mid) else (lo, mid)
    return lo


need = {}
for name, T in (('kappa_x2', T2), ('kappa_0.01', T01)):
    rec = {'a_target': f(T)}
    # (i) max decoder charge C (as rank-8 copies) at u=1 and u=2, everything else fixed
    for u in (1, 2):
        cmax = bisect(lambda c: ok(plus(data, helper, {8: Q(c) / 8} if c else {}), T, centre=Q(c), unpaid=u), Q(0), Q(400))
        rec[f'max_C_u{u}'] = None if cmax is None else f(cmax, 5)
    # (ii) minimal number of unpaid coordinates u (fractional) at C=64 and C=0
    for cc in (64, 0):
        loc = LOCAL if cc else plus(data, helper)
        umin = bisect(lambda s: not ok(loc, T, centre=cc, unpaid=1 + s), Q(0), Q(20))
        rec[f'min_u_at_C{cc}'] = f(1 + umin, 5) if umin is not None else '<=1'
    # (iii) max surviving fraction of helper climbs (same rank mix), rest idle, u=1 C=64
    h = bisect(lambda t: ok(plus(data, centre, {r: n * t for r, n in helper.items()}), T), Q(0), Q(1))
    rec['max_helper_fraction_u1'] = None if h is None else f(h, 4)
    # (iv) minimal uniform child rank for all non-centre local mass (u=1, C=64)
    lm = sum(r * n for r, n in plus(data, helper).items())
    rr = next((r for r in range(1, 52) if ok(plus({r: Q(lm, r)}, centre), T)), None)
    rec['min_uniform_local_rank_u1'] = rr
    rr2 = next((r for r in range(1, 52) if ok(plus({r: Q(lm, r)}), T, centre=0)), None)
    rec['min_uniform_local_rank_u1_C0'] = rr2
    need[name] = rec
out['requirements'] = need

# combined levers: max decoder charge C (as rank-8 copies) on top of a denominator change
combo = {}
for lname, base_local in (('frame one-shot (data+helpers)', plus(data1, helper1)),
                          ('helper climbs removed', dict(data)),
                          ('helper climbs removed + data one-shot', dict(data1)),
                          ('width-9 helpers removed', plus(data, {r: n for r, n in helper.items()}))):
    if lname == 'width-9 helpers removed':
        # rebuild helper histogram without the 240 width-9 helpers
        hk = Counter()
        for ev in prog['forward']:
            if ev[0] == 'climb' and ev[2] >= 224 and widths_by_reg[ev[2]] != 9:
                hk[ev[5]] += 1
        base_local = plus(data, hk)
    rec = {'a_cap_at_C64_u1': f(cap(plus(base_local, centre)))}
    for name, T in (('kappa_x2', T2), ('kappa_0.01', T01)):
        for u in (1, 2):
            cmax = bisect(lambda c: ok(plus(base_local, {8: Q(c) / 8} if c else {}), T, centre=Q(c), unpaid=u), Q(0), Q(600))
            rec[f'{name}_max_C_u{u}'] = None if cmax is None else f(cmax, 5)
    combo[lname] = rec
out['combined'] = combo

# ---------------------------------------------------------------- 6. universal cap (any word, coverage u, halving guard R=51)
R = C.max_idle_rank(M73)
lnR = L(R)
for u in (1, 2):
    out[f'universal_cap_u{u}'] = f(Q(u) / ((M73 - u) * lnR[0]))

print(json.dumps(out, indent=1))
print('PASS_DEFICIT_LOGMASS_LAW')
