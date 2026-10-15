#!/usr/bin/env python3
"""Framework ceiling for the recursive-supplier family (planning arithmetic, not a construction).

Master formula (Jensen form of the tilt identity, an UPPER bound on the root a*):

    a*  <=  A := ln( m / (m - nu) ) / lbar,          kappa <= A / (1 + A)
    nu   = u - s/(s-1) * C / V                        (net unpaid coordinates per family)
    m    = u + s (h - 1)                              (template identity; selected 73 = 1 + 9*8)
    lbar = phi * ln(m / r_loc) + (1 - phi) * ln(m / R),   R = max rank with m^2 > 2 R^2

phi = share of paid mass that the local word forces at rank <= r_loc (data/helper climbs),
the rest is idle mass packable up to the halving-guard rank R.

Derivation: Psi(a) = W^-1 sum n_r (r/m)^(1-a) = 1  <=>  E_pi exp(a L) = 1 + D/M  with
L_r = ln(m/r), pi_r = r n_r / M (tilt identity).  Jensen: exp(a lbar) <= 1 + D/M, and per family
D/W = nu, M/W = m - nu, so a* <= ln(m/(m-nu))/lbar.  lbar is minimised by putting each class at
its largest admissible rank.  Everything else (helpers, frames, hosting) only raises lbar or
lowers nu, so A is a ceiling for every construction obeying the stated identities.

usage: python3 -I -B framework_ceiling.py
Standard library only.
"""
import math

KAPPA_NOW = 0.0017            # our PR #2679 ~0.0017149; selected 0.0016727
KAPPA_SEL = 0.001672660197076918


def R_guard(m):
    R = int(m / math.sqrt(2))
    while 2 * (R + 1) ** 2 < m * m:
        R += 1
    while 2 * R * R >= m * m:
        R -= 1
    return max(R, 1)


def nu_of(u, s, C, V):
    return u - s / (s - 1) * C / V


def A_cap(m, nu, lbar):
    if nu <= 0:
        return 0.0
    return math.log(m / (m - nu)) / lbar


def kap(a):
    return a / (1 + a)


def lbar_of(m, phi, r_loc, R=None):
    R = R_guard(m) if R is None else R
    return phi * math.log(m / r_loc) + (1 - phi) * math.log(m / R)


def regimes(m, s, h, phi_word=0.7724, r_word=1.899):
    """lbar for three rank regimes: word-like (selected local word), data-floor (only the forced
    data climbs, share 2/(s-1), at one-shot rank h-1; everything else ideal at R), coverage-only (all at R)."""
    pd = min(1.0, 2 / (s - 1))
    return [lbar_of(m, phi_word, r_word), lbar_of(m, pd, h - 1), lbar_of(m, 0.0, h)]


def exact_root(hist, m, W):
    """Root of sum n_r (r/m)^(1-a) = W by bisection (float)."""
    f = lambda a: sum(n * (r / m) ** (1 - a) for r, n in hist.items()) - W
    lo, hi = 0.0, 1.0
    assert f(lo) < 0 < f(hi)
    for _ in range(200):
        mid = (lo + hi) / 2
        (lo, hi) = (mid, hi) if f(mid) < 0 else (lo, mid)
    return lo


def fmt(a):
    k = kap(a)
    return f"a<={a:.5f} kappa<={k:.5f} (x{k / KAPPA_NOW:6.1f})"


def section(t):
    print("\n" + "=" * 100 + "\n" + t + "\n" + "=" * 100)


def main():
    # ------------------------------------------------------------------ A. calibration
    section("A. Calibration on the selected record (#1094 complex side)")
    G = {1: 349604, 2: 60192, 3: 62585, 4: 20448, 5: 11808, 7: 7, 8: 1152, 48: 4866}
    m, W = 73, 14336
    a_true = exact_root(G, m, W)
    M = sum(r * n for r, n in G.items())
    D = m * W - M
    lb = sum(r * n * math.log(m / r) for r, n in G.items()) / M
    print(f"exact root a*={a_true:.10f}  (README: 0.001675462677)")
    print(f"nu=D/W={D / W:.6f} (5/14={5 / 14:.6f})  lbar={lb:.4f}  r_eff={m * math.exp(-lb):.3f}")
    print(f"Jensen cap ln(1+D/M)/lbar = {math.log(1 + D / M) / lb:.10f}   nu/Lambda = {D / W / (M / W * lb):.10f}")
    # per-call local word + ideal idle packing (complex-local-word-ceiling)
    word = {1: 2412, 2: 418, 3: 432, 4: 142, 5: 82, 8: 8}
    w = 896 / 9
    loc = sum(r * n for r, n in word.items())
    tot = 72 * w + 64
    phi_sel = loc / tot
    lb_loc = sum(r * n * math.log(73 / r) for r, n in word.items()) / loc
    r_loc_sel = 73 * math.exp(-lb_loc)
    hist = dict(word)
    hist[51] = (tot - loc) / 51
    print(f"selected word ceiling (idle at R=51), exact root: a={exact_root(hist, 73, w):.10f}  (README 0.001690054)")
    print(f"local share phi={phi_sel:.4f}, local effective rank r_loc={r_loc_sel:.3f}")
    nu_sel = nu_of(1, 9, 64, 112)
    print("master formula, selected parameters:", fmt(A_cap(73, nu_sel, lbar_of(73, phi_sel, r_loc_sel))))

    # ------------------------------------------------------------------ B. levers inside the template
    section("B. One- and multi-lever ceilings inside the selected template (m=73, s=9, h=9, V=112)")
    rows = [
        ("selected (u=1, C=64, local word as-is)", 1, 64, phi_sel, r_loc_sel),
        ("local word all at rank h=9", 1, 64, phi_sel, 9),
        ("helpers deleted: data climbs only (phi=.25) at rank 8", 1, 64, 0.25, 8),
        ("no local word at all (phi=0, all mass at R=51)", 1, 64, 0.0, 9),
        ("C=0 (decoder free), local as-is", 1, 0, phi_sel, r_loc_sel),
        ("C=0 and phi=0 (every lever free, u=1)", 1, 0, 0.0, 9),
        ("u=2, C=128 (plane ports at line ratio)", 2, 128, phi_sel, r_loc_sel),
        ("u=2, C=224 (best known plane decoder)", 2, 224, phi_sel, r_loc_sel),
        ("u=4=n/2 (orthogonality max), C=0, phi=0", 4, 0, 0.0, 9),
    ]
    for name, u, C, phi, rl in rows:
        mm = 73 + (u - 1)  # m = u + s(h-1)
        nu = nu_of(u, 9, C, 112)
        print(f"{name:58s} nu={nu:6.3f} m={mm:3d} lbar={lbar_of(mm, phi, rl):.3f}  {fmt(A_cap(mm, nu, lbar_of(mm, phi, rl)))}")

    # ------------------------------------------------------------------ C. stage count s
    section("C. Stage count s (E8 lines u=1, h=9, V=112, C=64; m=1+8s) -- if s were a free template parameter")
    print("data floor: 2V data registers per call each climb h-1 inside the stage => phi_data = 2/(s-1) (=0.25 at s=9, exact)")
    for s in [3, 4, 5, 6, 7, 9, 12, 16]:
        mm = 1 + 8 * s
        nu = nu_of(1, s, 64, 112)
        cells = [kap(A_cap(mm, nu, lb_)) for lb_ in regimes(mm, s, 9)]
        print(f"s={s:2d} m={mm:3d} nu={nu:6.3f} eps=nu/m={nu / mm:.5f}  kappa caps [word | data-floor | coverage-only]: "
              + " ".join(f"{c:.5f}" for c in cells))

    # ------------------------------------------------------------------ D. lattice / code families, u=1
    section("D. Line-port codes, u=1 (V lines, C floor), h=n+1; s=n+1 and best s; regimes word / data-floor / coverage-only")
    codes = [  # name, n, V (lines), C floor used
        ("D4", 4, 12, 12), ("E6", 6, 36, 36), ("E7 (gain 14)", 7, 63, 49), ("E8 P112", 8, 112, 64),
        ("E8 P120 (gain 48)", 8, 120, 72), ("E8 P112 C=56 (open lever)", 8, 112, 56),
        ("BW16 (C>=n^2, heuristic)", 16, 2160, 256), ("Leech (C>=n^2, heuristic)", 24, 98280, 576),
        ("C=0 (no decoder cost) n=8", 8, 112, 0),
    ]
    for name, n, V, C in codes:
        h = n + 1
        def caps(s):
            mm = 1 + s * n
            nu = nu_of(1, s, C, V)
            return [kap(A_cap(mm, nu, lb_)) for lb_ in regimes(mm, s, h)]
        c1 = caps(n + 1)
        best = max(range(3, 80), key=lambda s: caps(s)[1])
        c2 = caps(best)
        print(f"{name:28s} s=n+1: " + " / ".join(f"{c:.5f}" for c in c1)
              + f"   | best s={best:2d}: " + " / ".join(f"{c:.5f}" for c in c2))

    # ------------------------------------------------------------------ E. growth directions
    section("E. Growth directions: h -> infinity (data-floor regime unless stated)")
    print("E1. u = 1, gamma = C/V = 4/7 fixed, s=9:   kappa cap  [data-floor | coverage-only]")
    for h in [5, 9, 17, 33, 65, 129, 1025]:
        mm = 1 + 9 * (h - 1)
        nu = nu_of(1, 9, 4 / 7, 1)
        r = regimes(mm, 9, h)
        print(f"  h={h:5d} m={mm:6d}  {kap(A_cap(mm, nu, r[1])):.6f}  {kap(A_cap(mm, nu, r[2])):.6f}")
    print("E2. port labels of dimension u = k(h-1) (k<=1/2: orthogonal neighbours need dims summing <= h-1), m = u + s(h-1)")
    print("    gamma = 0 (free decoder) / 1/2 (Haemers floor when k=1/2) / 4/7 (E8 line ratio); kappa caps, data-floor regime")
    for k in [0.5, 0.25]:
        for s in [3, 4, 5, 9]:
            for h in [9, 129, 10 ** 6 + 1]:
                u = k * (h - 1)
                mm = u + s * (h - 1)
                cells = []
                for g in [0, 0.5, 4 / 7]:
                    nu = u * (1 - s / (s - 1) * g)
                    cells.append(kap(A_cap(mm, nu, regimes(mm, s, h)[1])))
                print(f"  k={k} s={s} h={h:8d} eps(g=0)={u / mm:.4f}  " + " / ".join(f"{c:.4f}" for c in cells))

    # ------------------------------------------------------------------ G. global maximum
    section("G. Global maximum of the data-floor ceiling over (s, h, u<=h/2, gamma)")
    for g, lbl in [(0.0, "gamma=0 (decoder free)"), (0.5, "gamma=1/2"), (4 / 7, "gamma=4/7 (E8 line ratio)")]:
      for hmin in (3, 9):
        best = None
        for s in range(3, 40):
            for h in range(hmin, 200):
                for u in range(1, h // 2 + 1):   # I6: pairwise-orthogonal u-spaces in an h-space (ctoth: <=4 planes in F2^9)
                    mm = u + s * (h - 1)
                    nu = u * (1 - s / (s - 1) * g)
                    a = A_cap(mm, nu, regimes(mm, s, h)[1])
                    if best is None or a > best[0]:
                        best = (a, s, h, u, mm, nu)
        a, s, h, u, mm, nu = best
        print(f"{lbl:28s} h>={hmin}: max at s={s} h={h} u={u} m={mm} nu={nu:.3f}: {fmt(a)}")
    print("(u=1 restricted, gamma=4/7, any s,h:)", end=" ")
    best = max(((A_cap(1 + s * (h - 1), nu_of(1, s, 4 / 7, 1), regimes(1 + s * (h - 1), s, h)[1]), s, h)
                for s in range(3, 40) for h in range(3, 200)))
    print(f"s={best[1]} h={best[2]}: {fmt(best[0])}")

    # ------------------------------------------------------------------ F. what epsilon a target needs
    section("F. Net unpaid fraction eps = nu/m needed for a target kappa  (eps >= 1 - exp(-a lbar))")
    for kt in [0.0034, 0.01, 0.1, 0.5]:
        a = kt / (1 - kt)
        cells = []
        for lbl, l in [("lbar=ln sqrt2 (all at R)", math.log(math.sqrt(2))), ("lbar=ln 3", math.log(3)),
                       ("lbar=2.09 (rank 9 at m=73)", math.log(73 / 9)), ("lbar=2.93 (today)", 2.9254)]:
            cells.append(f"{lbl}: {1 - math.exp(-a * l):.4f}")
        print(f"kappa={kt:<7} a={a:.4f}  " + " | ".join(cells))
    print(f"today eps = (5/14)/73 = {5 / 14 / 73:.5f}")
    print("kappa=1 is unreachable for every finite a (kappa=a/(1+a)<1); kappa->1 needs a->inf, i.e. eps->1.")

    # ------------------------------------------------------------------ H. decoder compression needed
    section("H. Largest decoder ratio gamma = C/(uV) allowed for a target, vs port-label fraction k = u/(h-1) (data-floor, h large)")
    h = 10 ** 4 + 1
    for kt in [0.0034, 0.01, 0.05, 0.1]:
        at = kt / (1 - kt)
        cells = []
        for k in [1 / 8, 1 / 4, 1 / 3, 0.45, 0.5]:
            best = None
            for s in range(3, 40):
                u = k * (h - 1)
                mm = u + s * (h - 1)
                lb_ = regimes(mm, s, h)[1]
                nu_need = mm * (1 - math.exp(-at * lb_))      # A >= at  <=>  nu >= m(1-exp(-at lbar))
                g = (1 - nu_need / u) * (s - 1) / s
                if best is None or g > best[0]:
                    best = (g, s)
            cells.append(f"k={k:.3f}: g<={best[0]:+.3f} (s={best[1]})")
        print(f"kappa>={kt:<6}  " + " | ".join(cells))
    print("closed-form supremum (gamma=0, k=1/2, s=3): A = ln(7/6)/ln(7/2) =", f"{math.log(7 / 6) / math.log(3.5):.5f}",
          "-> kappa <", f"{kap(math.log(7 / 6) / math.log(3.5)):.5f}")



if __name__ == "__main__":
    main()
