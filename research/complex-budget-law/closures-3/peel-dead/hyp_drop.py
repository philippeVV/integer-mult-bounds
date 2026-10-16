#!/usr/bin/env python3
"""expQH: which hypothesis of the framework ceiling must break for kappa -> 1?

Standard library only:  python3 -I -B hyp_drop.py   (a few seconds)

Part A  re-derives the family supremum of framework-ceiling/ and then drops ONE hypothesis at a time
        (I4 data floor, I5 halving guard, I6 u<=h/2, I7 s>=3, I3 template, transfer), keeping the rest.
Part B  the "peel" limit (top-ranked violation): children of rank close to m.  Exact tilt roots of
        explicit per-family histograms versus the closed forms
            a_peel  -> nu                                   (one climb per family, m -> infinity)
            a_touch <= ln(m/(m-nu)) * (m-u) / [ (c-1) ln m + (m-u-c+1) ln(m/(m-u-c+1)) ]
        (c touch-separated climbs per data family, Lemma 3 extreme point = most favourable).
Part C  the first precise obstruction: live-family factor (s-1)/s of the coverage law + c >= stages
        touches + decoder floor; best a over m for each (u=1, gamma, s).
Part D  width chain m -> m-1 -> ... (fixed m cannot peel; a chain gives the mediant, = nu in the limit).
Part E  outer caps: unfused a/(1+a), Gaussian range kappa < 1/4, fused ~ a.
"""
import math

K0 = 0.0017


def kap(a):
    return a / (1 + a)


def R_guard(m):
    R = int(m / math.sqrt(2))
    while 2 * (R + 1) ** 2 < m * m:
        R += 1
    while 2 * R * R >= m * m:
        R -= 1
    return max(R, 1)


def A_cap(m, nu, lbar):
    if nu <= 0 or lbar <= 0:
        return 0.0
    return math.log(m / (m - nu)) / lbar


def lbar_floor(m, s, h, data_floor=True, guard=True):
    R = R_guard(m) if guard else m - 1
    pd = min(1.0, 2 / (s - 1)) if data_floor else 0.0
    return pd * math.log(m / (h - 1)) + (1 - pd) * math.log(m / R)


def root(hist, m, W):
    """exact root of sum n_r (r/m)^(1-a) = W (bisection; a may exceed 1)."""
    f = lambda a: sum(n * (r / m) ** (1 - a) for r, n in hist.items()) - W
    lo, hi = 0.0, 1.0
    if f(lo) >= 0:
        return 0.0
    while f(hi) < 0:
        hi *= 2
        if hi > 1e6:
            return float('inf')
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
    return lo


def sec(t):
    print("\n" + "=" * 96 + "\n" + t + "\n" + "=" * 96)


def sup(gamma=0.0, hmin=3, smin=3, umax=lambda h: h // 2, template=True, data_floor=True, guard=True,
        hmax=200, smax=40):
    best = None
    for s in range(smin, smax):
        for h in range(hmin, hmax):
            for u in range(1, umax(h) + 1):
                m = u + s * (h - 1)
                nu = u * (1 - s / (s - 1) * gamma)
                a = A_cap(m, nu, lbar_floor(m, s, h, data_floor, guard))
                if best is None or a > best[0]:
                    best = (a, s, h, u, m)
    return best


def show(lbl, b):
    a0, s, h, u, m = b
    a = min(a0, 1.0)   # reading cost: cost/volume >= const, so a <= 1
    print(f"{lbl:58s} s={s:2d} h={h:3d} u={u:4d} m={m:5d}  A={a0:8.4f} a<={a:.4f} kappa<={kap(a):.4f} (x{kap(a)/K0:6.1f})"
          f"  fused~a (x{a/K0:6.1f})")


def main():
    sec("A. Drop one hypothesis at a time (free decoder gamma=0 unless stated); grid s<40, h<200")
    show("all kept (I3-I7), h>=3  [framework-ceiling: 0.134]", sup())
    show("all kept, h>=9          [framework-ceiling: 0.118]", sup(hmin=9))
    show("all kept, gamma=1/2 (Haemers floor)", sup(gamma=0.5))
    show("drop I5 halving guard (R=m-1)", sup(guard=False))
    show("drop I4 data floor (idle at R)", sup(data_floor=False))
    show("drop I7 (s=2 allowed)", sup(smin=2))
    show("drop I6 -> u<=h-1", sup(umax=lambda h: h - 1))
    show("drop I6 -> u<=h", sup(umax=lambda h: h))
    show("drop I6 entirely, u<=20h (labels outside stage space)", sup(umax=lambda h: 20 * h, hmax=40))
    show("drop I6 entirely, gamma=1/2", sup(gamma=0.5, umax=lambda h: 20 * h, hmax=40))
    b = sup(data_floor=False, guard=False)
    show("drop I4+I5 (all paid mass at rank m-1)", b)
    print("   NB: a>1 means fewer than one child per role: cost/volume bounded => a is capped at 1 by reading cost.")
    print("   closed form u=h/2 s=3, h->inf, I6 kept: A=ln(7/6)/ln(7/2) =", round(math.log(7 / 6) / math.log(3.5), 5))
    for x in [0.5, 1, 2, 5, 10, 100, 1e4]:
        A = math.log((x + 3) / 3) / math.log(x + 3)
        print(f"   I6 dropped, s=3, x=u/(h-1)={x:<7}: A=ln((x+3)/3)/ln(x+3)={A:.4f} kappa<={kap(A):.4f} (x{kap(A)/K0:.0f})")

    sec("B. Peel limit: per family (c-1) unit climbs + one big climb (Lemma 3 extreme) + decoder mass u-nu")
    print(" exact tilt roots, W=1 family, fractional counts. dec@big = decoder children at the big rank (favourable);")
    print(" dec@1 = decoder children at rank 1 (pessimistic). Closed form: a ~ nu/(1 + mu ln m), mu = small-rank mass/family")
    for nu, u in [(5 / 14, 1), (1 / 7, 1), (0.5, 1), (1.0, 1)]:
        for c in [1, 2, 3, 9]:
            row = []
            for m in [9, 73, 1000, 10 ** 5]:
                big = m - u - (c - 1)
                if big < 1:
                    row.append("        -        ")
                    continue
                out = []
                for rdec in (big, 1):
                    hist = {}
                    hist[1] = hist.get(1, 0) + (c - 1)
                    hist[big] = hist.get(big, 0) + 1
                    if u - nu > 0:
                        hist[rdec] = hist.get(rdec, 0) + (u - nu) / rdec
                    out.append(root(hist, m, 1))
                mu = (c - 1) + ((u - nu) if True else 0)
                row.append(f"{out[0]:.4f}|{out[1]:.4f}")
            print(f" nu={nu:.4f} c={c}: m=9,73,1e3,1e5 -> dec@big|dec@1: " + "  ".join(row))
    print(" reading: one climb per family and decoder at large rank -> a -> nu (m-independent);")
    print("          each unit of per-family mass at O(1) rank costs ~ln m: a ~ nu/(1 + mu ln m).")
    for nu, mu in [(5 / 14, 0.0), (5 / 14, 1.0), (5 / 14, 2.0), (1.0, 2.0)]:
        print(f"   closed form nu={nu:.3f} mu={mu}: m=73 -> {nu/(1+mu*math.log(73)):.4f}, m=1e5 -> {nu/(1+mu*math.log(1e5)):.4f}")


    sec("B2. Lemma P check: a <= ln[(m/B)(1 + (nu + (u-nu)(u+mu)/B)/mu)] / ln m,  B = m-u-mu (decoder at rank B)")
    worst = 0.0
    n = 0
    for u in [1, 2, 4]:
        for nu in [0.1, 5 / 14, 0.5, 1.0]:
            if nu > u:
                continue
            for mu in [1, 2, 3, 8]:
                for m in [u + mu + 2, 12, 73, 500, 10 ** 4, 10 ** 6]:
                    B = m - u - mu
                    if B < 1:
                        continue
                    hist = {1: mu}
                    hist[B] = hist.get(B, 0) + 1 + (u - nu) / B
                    a = root(hist, m, 1)
                    bound = math.log((m / B) * (1 + (nu + (u - nu) * (u + mu) / B) / mu)) / math.log(m)
                    assert a <= bound + 1e-12, (u, nu, mu, m, a, bound)
                    worst = max(worst, a / bound)
                    n += 1
    print(f" {n} cases: exact root <= Lemma P bound in all; max root/bound = {worst:.4f}")
    for mu in [1, 2]:
        print(f" nu=5/14 u=1 mu={mu}: m=73 bound {math.log((73/(72-mu))*(1+(5/14+(9/14)*(1+mu)/(72-mu))/mu))/math.log(73):.4f}"
              f", m=1e6 bound {math.log((1e6/(1e6-1-mu))*(1+(5/14)/mu))/math.log(1e6):.4f}")

    sec("C. First obstruction: coverage live-family factor + touch count c >= s (data touched every stage)")
    print(" nu = u (1 - s/(s-1) gamma) [I2];  data touched >= s times => c >= s;  best m (u=1); decoder mass charged at")
    print(" the climbs' mean log-rank (favourable Jensen form)")
    for gname, g in [("gamma=0", 0.0), ("gamma=7/13 (8-term floor)", 7 / 13), ("gamma=4/7 (E8 lines)", 4 / 7)]:
        for s in [3, 4, 5, 9]:
            nu = 1 - s / (s - 1) * g
            best = (0, None)
            for m in range(s + 1, 5000):
                big = m - 1 - (s - 1)
                lb = ((s - 1) * math.log(m) + big * math.log(m / big)) / (m - 1)
                A = A_cap(m, nu, lb)
                if A > best[0]:
                    best = (A, m)
            A, m = best
            print(f" {gname:26s} s={s}: nu={nu:+.4f}  best m={m}  a<={A:.4f} kappa<={kap(A):.4f} (x{kap(A)/K0:5.1f})")
    print(" => with u=1 the peel is pushed back to tiny m (small-unit lane NF/NG); at m=73,s=9 it is:")
    m, s, nu = 73, 9, 5 / 14
    big = m - 1 - (s - 1)
    lb = ((s - 1) * math.log(m) + big * math.log(m / big)) / (m - 1)
    A = A_cap(m, nu, lb)
    print(f"    lopsided data (8 unit + 1x{big}) helpers gone: a<={A:.4f} kappa<={kap(A):.4f} (x{kap(A)/K0:.1f})")
    print(" Lemma W (expQC): every restored helper/fixed register climbs its full m => eps <= 1 - |Fix|/W.")

    sec("D. Width chain: fixed m cannot peel (r=m-1 is coprime to m); chain m_k = k, nu per level")
    for nu in [5 / 14, 1.0]:
        for K in [10, 100, 10 ** 4]:
            num = sum(math.log(k / (k - nu)) for k in range(2, K + 1) if k > nu)
            den = sum(math.log(k / (k - 1)) for k in range(2, K + 1))
            print(f" nu={nu:.4f} widths 2..{K}: mediant a = {num / den:.4f}   (-> nu = {nu:.4f})")

    sec("E. Outer caps (assembly-domain-audit.md:22-26)")
    for a in [0.0016755, 0.123, 0.357, 1.0]:
        print(f" a={a:<9} unfused kappa<=a/(1+a)={kap(a):.4f}; Gaussian range kappa<1/4 -> {min(kap(a), .25):.4f}"
              f"; fused ~a -> {min(a, 1):.4f} (x{min(a, 1)/K0:.0f})")


if __name__ == "__main__":
    main()
    print("PASS_EXPQH_HYPDROP")  # every Lemma P case is asserted inside main()
