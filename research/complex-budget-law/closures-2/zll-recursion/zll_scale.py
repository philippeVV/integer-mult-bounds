#!/usr/bin/env python3
"""expQC: does a ZLL-style recursive construction help as the bit supplier's local word?

Stdlib only, < 1 s.  Four checks:
  A. Tensor invariance of the tilt-identity root a (supplier recursion == tensor power),
     plus the bracket a(A(x)B) in [min, max] for level-dependent seeds (the mediant formula is exact only for
     single-rank levels; see check_qg.py Part 2 and README).
  B. HWZ auxiliary recursion counts used by ZLL (h = 0..3, u = 58): sessions per edge k/m,
     and the two ceilings this gives for ZLL's own compiled circuit:
       Lemma W (fixed-point dilution):   delta <= 2k/N <= 2k/(3k+2m)
       ZLL temporal certificate scale:   eps0 ~ 2k/(N(M+1)) <= 2k/((3k+2m)*12m)
  C. ZLL tensor family S_j compiled once at depth j: per-seed exponent is constant
     (ln(1/rho_j)/ln b_j = ln(1/rho_0)/ln b_0), while the compiled deficit ceiling 2k_j/(3k_j+2m_j)
     decays doubly exponentially because m_j/k_j = (m_0/k_0)^(2^j).
  D. Uniform-metric dictionary a = 1 - ln(m/k)/ln b on known seeds.
"""
import math
from fractions import Fraction as Fr

K0 = 0.0017  # reference kappa for multipliers


def root_a(hist, m, W, lo=0.0, hi=0.999):
    """Root of sum n_r (r/m)^(1-a) = W (bisection)."""
    f = lambda a: sum(n * (r / m) ** (1 - a) for r, n in hist.items()) - W
    if f(lo) >= 0:
        return 0.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return lo


def tensor(h1, m1, W1, h2, m2, W2):
    h = {}
    for r1, n1 in h1.items():
        for r2, n2 in h2.items():
            h[r1 * r2] = h.get(r1 * r2, 0) + n1 * n2
    return h, m1 * m2, W1 * W2


def part_A():
    print("== A. tensor invariance of the supplier root (tilt identity) ==")
    # selected complex global profile (complex-budget-law README); m solved from a* = 0.0016754627
    prof = {1: 349604, 2: 60192, 3: 62585, 4: 20448, 5: 11808, 7: 7, 8: 1152, 48: 4866}
    W = 14336
    a_star = 0.001675462676808
    m = (sum(n * r ** (1 - a_star) for r, n in prof.items()) / W) ** (1 / (1 - a_star))
    a1 = root_a(prof, m, W)
    h2, m2, W2 = tensor(prof, m, W, prof, m, W)
    a2 = root_a(h2, m2, W2)
    print(f"  selected profile: implied m = {m:.4f}, root a = {a1:.12f}")
    print(f"  self-tensor (two supplier levels as ONE seed): root a = {a2:.12f}  (diff {a2-a1:+.1e})")
    # synthetic seeds with different exponents -> tensor root lies between them (not a mediant in general)
    sA = ({1: 30, 2: 10}, 8, 9)
    sB = ({3: 5, 1: 20}, 6, 6)
    aA, aB = root_a(*sA), root_a(*sB)
    hAB, mAB, WAB = tensor(*sA, *sB)
    aAB = root_a(hAB, mAB, WAB)
    print(f"  seed A a = {aA:.6f}, seed B a = {aB:.6f}, A(x)B a = {aAB:.6f}  "
          f"(between: {min(aA,aB) <= aAB <= max(aA,aB)})")
    return a1, a2


def hwz_level(h, u):
    """Return (R=k/f, U=m/f, a=duration (full duplex), b=deleted distance) for HWZ level h, target u.
    ZLL Appendix B: base q=ceil(u): k=q, f=1, m=3q^2-q+1, a=3, b=5.
    Level h+1 at target u: outer = level h at target 3u, inner = level h at target max(5, U1).
      d = 2(m1-f1), nR/nL = d/k2;  k+ = nL k1, f+ = nL f1 + nR f2, m+ = a2 nL f1 + nR m2.
    Ratios only (cover sizes cancel)."""
    if h == 0:
        q = math.ceil(u)
        return Fr(q), Fr(3 * q * q - q + 1), 3, 5
    R1, U1, a1, b1 = hwz_level(h - 1, 3 * u)
    R2, U2, a2, b2 = hwz_level(h - 1, max(Fr(5), U1))
    # per unit f1 nL:  k = R1, f = 1 + (nR/nL) f2/f1, m = a2 + (nR/nL) m2/f1
    # nR/nL = d/k2 = 2 f1 (U1 - 1)/k2  ->  (nR/nL)/f1 * f2 = 2(U1-1)/R2 ; (nR/nL)/f1 * m2 = 2(U1-1) U2/R2
    t = 2 * (U1 - 1) / R2
    k = R1
    f = 1 + t
    m = a2 + t * U2
    return k / f, m / f, a1 * a2, b1 * b2


def part_B():
    print("\n== B. ZLL's own compiled circuit (HWZ auxiliary at h, u = 58) ==")
    rows = []
    for h in range(0, 4):
        R, U, a, b = hwz_level(h, Fr(58))
        mk = U / R  # edges per session
        lmk = math.log(mk.numerator) - math.log(mk.denominator)   # ln(m/k), big-int safe
        # Lemma W: delta <= 2/(3+2 m/k) ~ k/m ; temporal scale eps0 <= 2/((3+2m/k) 12 m/k) ~ (k/m)^2/12
        l10W = (math.log(2) - math.log(3 + 2 * math.exp(lmk)) if lmk < 700 else math.log(1) - lmk) / math.log(10)
        l10e = (l10W * math.log(10) - math.log(12) - lmk) / math.log(10)
        rows.append((h, float(R), lmk, a, b, l10W, l10e))
        print(f"  h={h}: k/f={float(R):.3e}  ln(m/k)={lmk:.4e}  T(full-duplex)=3^(2^{h})  b=5^(2^{h})  "
              f"log10 delta_W<={l10W:.4g}  log10 eps0<={l10e:.4g}")
    return rows


def part_C(eps0, b0, mk0, K=3):
    print(f"\n== C. tensor family S_j (ZLL Sec. V), seed eps0={eps0:.3g}, b0={b0:.3g}, m0/k0={mk0:.3g} ==")
    rho0 = 1 - eps0
    for j in range(0, 6):
        lg = (2 ** j) * math.log(1 / rho0)          # ln(1/rho_j)
        L = (2 ** j) * math.log(b0)                 # hop metric: x_e = 1/b_j -> L = ln b_j
        a_net = lg / L                              # if frames realised the network metric directly
        ln_mk = (2 ** j) * math.log(mk0)            # ln(m_j/k_j)
        dW = 2 / (3 + 2 * math.exp(min(ln_mk, 700)))  # compiled-once ceiling (Lemma W)
        dC = "dead (K*rho>=1)" if K * math.exp(-lg) >= 1 else "Lemma C passes"
        print(f"  j={j}: ln gap={lg:.3e}  L=ln b_j={L:.3e}  a(net metric)={a_net:.3e}  "
              f"compiled delta<= {dW:.3e}  [{dC}]")


def part_D():
    print("\n== D. uniform hop-metric dictionary a = 1 - ln(m/k)/ln b (net used directly as wires) ==")
    seeds = [
        ("BH PG(2,9) core (dist 5, 735/157)", 735 / 157, 5),
        ("PG(2,16)-M nondet (dist 5, (4368-772)/772 physical)", (4368 - 772) / 772, 5),
        ("HWZ base q=58 (completion-time gadget; dist 3 WITH cut edge)", (3 * 58 * 58 - 58 + 1) / 58, 3),
    ]
    for name, mk, b in seeds:
        a = 1 - math.log(mk) / math.log(b)
        print(f"  {name}: a = {a:+.4f}  -> {'no gap (rho>=1)' if a <= 0 else f'x{a/K0:.0f}'}")


if __name__ == "__main__":
    part_A()
    rows = part_B()
    # Seed for part C: ZLL's normalized seed from HWZ h=3: eps0 ~ certificate scale, b0 ~ 2(M+3)(m*+1)
    # Seed for part C: HYPOTHETICAL generous stand-in using the h=0 counts (h=0 does not meet ZLL margin 24T<b;
    # the real seed needs h=3, whose eps0 is ~1e-1936, see B).
    h, R, lmk, a, b, l10W, l10e = rows[0]
    eps = 10 ** l10e
    mk = math.exp(lmk)
    M_over_k = 12 * mk                      # cross gates >= 12 m  (per designated session)
    b0 = 2 * (5 * M_over_k) * M_over_k      # b0 = 2h* ~ 2(m*+1)(M+3), m* >= 5M ; lower bound (k>=1)
    part_C(eps, b0, b0 * (1 - eps))
    part_D()
    a_seed = eps / math.log(b0)
    print(f"\nZLL seed built on the h=0 base (integer; killed by stronger-rank-screens sec.1). Ignoring sec.1:")
    print(f"  a <= eps0/ln b0 = {a_seed:.3e}  -> x{a_seed/K0:.1e} over 0.0017;  ZLL's actual h=3 seed: "
          f"log10 eps0 <= {rows[-1][6]:.4g}")
    print("PASS_EXPQC")
