#!/usr/bin/env python3
"""expNB: end-to-end kappa for a PAIR of suppliers built from one cross-field unit.

kappa_pair(h, x_bit, x_cx, dh=1, model) = A/(1+A),  A = min(a_bit, b_cx)
  bit supplier    : window h_b = h,      copy fraction x_bit, m_b = 1 + s_b (h_b - 1)
  complex supplier: window h_c = h + dh, copy fraction x_cx,  m_c = 1 + s_c (h_c - 1)
  (E8: h_b = 8 = rank_Q G, h_c = 9 = rank_F2 (J+G); dh = +1 is the odd-form price of
   complex labels for an even (bit-gain-capable) Gram, see README section 3.)

Per side the saving is modelled from the PR595 format ceiling X(s, h, x) (format-limits A.2,
re-implemented in float as in $S/expM2/format_map.py; validated 6.0005e-3 / 5.1863e-3):
  'ceiling'  : a = X*(h, x) = max_s X(s, h, x)            (theorem-level cap, eta = 1)
  'cal'      : a = a_E8 * X*(h, x) / X*(h_E8, x_E8)       (today's realised/ceiling ratio, best s)
  'cal9'     : a = a_E8 * X(9, h, x) / X(9, h_E8, x_E8)   (same, s pinned at 9 as in the words)
  'nohelp'   : a = eta_C * X*(h, x), eta_C = format-limits A.5 step C / step A at s = 9
               (real data fragmentation, every helper free; an optimistic realisation)
Calibration: a_bit(E8) = 0.0017334270 (#589), b(E8) = 0.0017149/(1-0.0017149) (our #2679),
so kappa_pair(8, 63/111, 64/112, model='cal') = 0.0017149.
Stdlib only, < 2 s.  Prints PASS_KAPPA_PAIR at the end.
"""
import math, sys
from functools import lru_cache

REF = 0.0017
KAPPA_2679 = 0.0017149
A_BIT_E8 = 0.0017334270                       # selected bit endpoint (#589), h=8, x=63/111
B_CX_E8 = KAPPA_2679 / (1 - KAPPA_2679)        # complex endpoint of #2679, h=9, x=64/112
E8_BIT = (8, 63 / 111)
E8_CX = (9, 64 / 112)
ETA_C = {'bit': 3.7868 / 4.6072, 'cx': 3.1544 / 3.9946}   # format-limits A.5 step C / step A


# ---------------- PR595 skeleton (format-limits A.2), k = 1 ----------------
def idle_slots(K, a):
    sl = []
    if 2 * K - 2 > 0: sl.append(a * (2 * K - 2))
    for p in range(1, K):
        sl += [2 * p * a, a]
        if 2 * K - 2 * p - 2 > 0: sl.append(a * (2 * K - 2 * p - 2))
    return sl


def slots(s, h, x):
    """forced slots per port: list of (cap, mass); m; deficit D_s"""
    K = (s - 1) // 2; a = h - 1; m = 1 + s * a
    out = [(a, 2 * s * a), (h, s * x)]
    out += [(max(r, h), 2 * r) for r in idle_slots(K, a)]
    return [(c, w) for c, w in out if w > 0 and c < m], m, (s - 1) - s * x


@lru_cache(maxsize=None)
def X(s, h, x):
    """root q of sum mass*((m/cap)^q - 1) = D_s  (the format ceiling for one supplier)"""
    sl, m, D = slots(s, h, x)
    if D <= 0: return 0.0
    f = lambda q: sum(w * ((m / c) ** q - 1) for c, w in sl) - D
    lo, hi = 0.0, 1.0
    if f(hi) < 0: return 1.0
    for _ in range(64):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
    return lo


def X_lin(s, h, x):
    """Jensen/e^y>=1+y upper bound D_s / Cost_s (closed form)"""
    sl, m, D = slots(s, h, x)
    return max(D, 0) / sum(w * math.log(m / c) for c, w in sl)


@lru_cache(maxsize=None)
def Xstar(h, x):
    # X(s) is unimodal in s with optimum s <= 7 on this grid (A.3 tail bound); scan to 41
    return max((X(s, h, x), s) for s in range(3, 42, 2))


# ---------------- per-side savings and the pair ----------------
def side_saving(side, h, x, model='cal'):
    hE, xE = E8_BIT if side == 'bit' else E8_CX
    aE = A_BIT_E8 if side == 'bit' else B_CX_E8
    if model == 'ceiling': return Xstar(h, x)[0]
    if model == 'cal': return aE * Xstar(h, x)[0] / Xstar(hE, xE)[0]
    if model == 'cal9': return aE * X(9, h, x) / X(9, hE, xE)
    if model == 'nohelp': return ETA_C[side] * Xstar(h, x)[0]
    raise ValueError(model)


def kappa_pair(h, x_bit, x_cx, dh=1, model='cal'):
    """h = bit window; complex window h + dh.  Returns (kappa, A, binding side, a_bit, b_cx)."""
    a = side_saving('bit', h, x_bit, model)
    b = side_saving('cx', h + dh, x_cx, model)
    A = min(a, b)
    return A / (1 + A), A, ('bit' if a < b else 'cx'), a, b


def twin_from_complex(h_c, x_c, rule='equal'):
    """bit twin of a unit found in COMPLEX orientation (F2 labels of rank h_c, Q decoder):
    bit window = rank_Q of the decoder Gram >= h_c - 1 (E8: 8 = 9 - 1); copy fraction
    'equal' (E8: 72/120 on both sides) or 'star' (x_b = x_c * h_c/(h_c-1))."""
    x_b = x_c if rule == 'equal' else x_c * h_c / (h_c - 1)
    return h_c - 1, x_b


# ---------------- report ----------------
def main():
    ok = True
    print("== validation (format-limits A.3) ==")
    v1, v2 = X(5, 8, 63 / 111), X(5, 9, 64 / 112)
    print(f"X(5, bit 8, 63/111) = {v1:.4e} (6.0005e-3);  X(5, cx 9, 64/112) = {v2:.4e} (5.1863e-3)")
    ok &= abs(v1 - 6.0005e-3) < 2e-7 and abs(v2 - 5.1863e-3) < 2e-7
    print(f"X(9) bit = {X(9, 8, 63/111):.4e} (4.6072e-3), cx = {X(9, 9, 64/112):.4e} (3.9946e-3)")
    k, A, side, a, b = kappa_pair(8, 63 / 111, 64 / 112, 1, 'cal')
    print(f"calibration: kappa_pair(8, 63/111, 64/112) = {k:.7f} (target {KAPPA_2679}); "
          f"a_bit={a:.7f} b_cx={b:.7f}; binds: {side}")
    ok &= abs(k - KAPPA_2679) < 1e-9 and side == 'cx'
    for mdl in ('cal9', 'nohelp', 'ceiling'):
        k, A, side, a, b = kappa_pair(8, 63 / 111, 64 / 112, 1, mdl)
        print(f"   E8 under model {mdl:8s}: kappa = {k:.5f} (x{k/REF:5.2f}), binds {side}")
    print("realised/ceiling ratios: bit %.4f, cx %.4f (best s);  s=9: bit %.4f, cx %.4f" % (
        A_BIT_E8 / Xstar(*E8_BIT)[0], B_CX_E8 / Xstar(*E8_CX)[0],
        A_BIT_E8 / X(9, *E8_BIT), B_CX_E8 / X(9, *E8_CX)))

    print("\n== closed form at s = 3 (K = 1: no idle slots), one side ==")
    print("X_3(h,x) = root of 6a((m/a)^q-1) + 3x((m/h)^q-1) = 2-3x, a=h-1, m=3h-2;")
    print("  linear bound (2-3x)/(6a ln(m/a) + 3x ln(m/h));  fit 0.29(1-1.5x)/(h-1)")
    for h in (5, 6, 7, 8, 9):
        for x in (0.3, 0.45, 0.6):
            print(f"  h={h} x={x:.2f}: root {X(3,h,x):.5f}  lin {X_lin(3,h,x):.5f}  fit {0.29*(1-1.5*x)/(h-1):.5f}"
                  f"   best-s X*={Xstar(h,x)[0]:.5f} (s={Xstar(h,x)[1]})")

    xs = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60]
    def table(title, fn):
        print("\n== " + title + " ==")
        print("h  " + "".join(f"  x={x:.2f}            " for x in xs))
        for h in (5, 6, 7):
            row = []
            for x in xs:
                k, A, side, a, b = fn(h, x)
                row.append(f"{k:.5f} x{k/REF:5.2f} {side:3s}")
            print(f"{h}  " + "  ".join(row))

    for mdl in ('cal', 'cal9', 'nohelp', 'ceiling'):
        table(f"[{mdl}] E8-like pair: bit (h, x), complex (h+1, x)",
              lambda h, x, mdl=mdl: kappa_pair(h, x, x, 1, mdl))
    for mdl in ('cal', 'ceiling'):
        table(f"[{mdl}] unit found in COMPLEX orientation at window h: complex (h, x), bit twin (h-1, x)",
              lambda h, x, mdl=mdl: kappa_pair(h - 1, x, x, 1, mdl))
        table(f"[{mdl}] sensitivity: bit twin pays x*h/(h-1) (bit copy +1/(h-1) relative)",
              lambda h, x, mdl=mdl: kappa_pair(h - 1, x * h / (h - 1), x, 1, mdl))
        table(f"[{mdl}] dh = 0 (both windows h, same x; upstream-triple-like)",
              lambda h, x, mdl=mdl: kappa_pair(h, x, x, 0, mdl))

    print("\n== new unit on ONE side only (the other keeps E8) ==")
    for mdl in ('cal', 'ceiling'):
        for h in (5, 6, 7):
            x = 0.4
            a_new = side_saving('cx', h, x, mdl); a_old = side_saving('bit', 8, 63 / 111, mdl)
            A = min(a_new, a_old)
            print(f"  [{mdl:7s}] complex new (h={h}, x={x}) b={a_new:.5f}, bit E8 a={a_old:.5f} "
                  f"-> kappa {A/(1+A):.5f} (x{A/(1+A)/REF:.2f})")

    print("\n== balance: x_bit at which the bit side stops binding, complex (h+1, x_cx), model cal ==")
    for h in (5, 6, 7, 8):
        for xc in (0.3, 0.4, 0.5):
            b = side_saving('cx', h + 1, xc, 'cal')
            lo, hi = 0.0, 1.0
            for _ in range(50):
                mid = (lo + hi) / 2
                if side_saving('bit', h, mid, 'cal') > b: lo = mid
                else: hi = mid
            print(f"  h_b={h} h_c={h+1} x_cx={xc:.2f}: bit binds iff x_bit > {lo:.4f}"
                  f" (ratio {lo/xc:.3f})")

    print("\n== x thresholds for a target multiplier (pair, E8-like dh=1, x_bit = x_cx) ==")
    for mdl in ('cal', 'ceiling'):
        for h in (4, 5, 6, 7, 8):
            out = []
            for mult in (2, 3, 5, 10):
                if kappa_pair(h, 0, 0, 1, mdl)[0] < mult * REF: out.append(' none '); continue
                lo, hi = 0.0, 1.0
                for _ in range(50):
                    mid = (lo + hi) / 2
                    if kappa_pair(h, mid, mid, 1, mdl)[0] >= mult * REF: lo = mid
                    else: hi = mid
                out.append(f"{lo:.4f}")
            print(f"  [{mdl:7s}] h={h}: x for x2/x3/x5/x10 = {' '.join(out)};"
                  f" kappa(x=0) = {kappa_pair(h,0,0,1,mdl)[0]:.4f}")


    print("\n== contract data for one supplier (per B_K word, T replicas) ==")
    print("m = 1 + s(h-1); W = (s-1) v T families; D = T[(s-1)v - s*cst]; rank sum = mW - D; nu = D/W = 1 - s x/(s-1)")
    for name, h, v, cst, s, T in (("bit E8 #555", 8, 111, 63, 9, 32), ("cx E8 #568", 9, 112, 64, 9, 8),
                                  ("cx new h=5", 5, 50, 20, 3, 1), ("bit twin h=4", 4, 50, 20, 3, 1)):
        m = 1 + s * (h - 1); W = (s - 1) * v * T; D = T * ((s - 1) * v - s * cst)
        a = X(s, h, cst / v)
        print(f"  {name:13s}: h={h} v={v} cst={cst} x={cst/v:.4f} s={s}: m={m} W={W} D={D} rank sum={m*W-D}"
              f" nu={D/W:.4f}; s_rank/W={(m*W-D)/W:.3f}; upstream rank-1 saving 1-log_m(s/W)={1-math.log((m*W-D)/W, m):.5f};"
              f" pressure-identity ceiling X(s)={a:.5f}")

    print("\n== expN2 LP floors on x (complex orientation, standard sinks unless noted) through the pair ==")
    for h, x, note in ((3, 9/11, 'std exact'), (3, 5/7, 'complement exact'), (4, 0.7822, 'std exact'),
                       (4, 0.4840, 'complement, local search'), (5, 0.5485, 'std exact'), (6, 0.4862, 'std partial')):
        kc = kappa_pair(h - 1, x, x, 1, 'cal')[0]; kx = kappa_pair(h - 1, x, x, 1, 'ceiling')[0]
        print(f"  h_c={h} x>={x:.4f} ({note}): pair kappa <= cal {kc:.5f} (x{kc/REF:.2f}), ceiling {kx:.5f} (x{kx/REF:.2f})")

    print("\nPASS_KAPPA_PAIR" if ok else "\nFAIL_KAPPA_PAIR")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
