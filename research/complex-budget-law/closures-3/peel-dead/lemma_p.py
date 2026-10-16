#!/usr/bin/env python3
"""expQQ: independent check of Lemma P (expQH sec.4) and of its sharpened form.

Family histogram (W = 1): mu unit climbs, one climb of rank B = m-u-mu, decoder mass u-nu at rank B.
  Psi(a) = mu m^(a-1) + (1 + (u-nu)/B) (B/m)^(1-a) = 1.
expQH bound :  a <= ln[(m/B)(1 + (nu + (u-nu)(u+mu)/B)/mu)] / ln m
sharpened P':  a <= ln(1 + nu/mu) / ln(m/r0)  -- for ANY placement of the remaining mass m-nu-mu at ranks <= m
               and the mu mass at ranks <= r0 (here r0 = 1).
Roots: 60-digit bisection with mpmath (independent of hyp_drop.py's float bisection).
"""
import math
import random
from mpmath import mp, mpf, log, power

mp.dps = 60


def root(hist, m, W=1):
    m = mpf(m)
    f = lambda a: sum(mpf(n) * power(mpf(r) / m, 1 - a) for r, n in hist) - W
    lo, hi = mpf(0), mpf(1)
    if f(lo) >= 0:
        return mpf(0)
    while f(hi) < 0:
        hi *= 2
    for _ in range(220):
        mid = (lo + hi) / 2
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return lo


def main():
    n = 0; worst_qh = 0; worst_sh = 0; viol = 0; sh_tighter = 0
    for u in [1, 2, 4]:
        for nu in [mpf('0.1'), mpf(5) / 14, mpf('0.5'), mpf(1)]:
            if nu > u:
                continue
            for mu in [1, 2, 3, 8]:
                for m in [u + mu + 2, 12, 73, 500, 10 ** 4, 10 ** 6]:
                    B = m - u - mu
                    if B < 1:
                        continue
                    hist = [(1, mu), (B, 1 + (u - nu) / B)]
                    a = root(hist, m)
                    bq = log((mpf(m) / B) * (1 + (nu + (u - nu) * (u + mu) / B) / mu)) / log(m)
                    bs = log(1 + nu / mu) / log(m)
                    n += 1
                    if a > bq or a > bs:
                        viol += 1
                        print('VIOLATION', u, nu, mu, m, a, bq, bs)
                    worst_qh = max(worst_qh, a / bq); worst_sh = max(worst_sh, a / bs)
                    sh_tighter += bs <= bq
    print(f'{n} cases (expQH grid). violations: {viol}. max root/bound: expQH {float(worst_qh):.4f}, '
          f'sharpened {float(worst_sh):.4f}; sharpened <= expQH bound in {sh_tighter}/{n}')
    # placement-free check of P': random histograms, mu mass at ranks <= r0, rest anywhere <= m
    rng = random.Random(7)
    worst = 0; k = 0
    for _ in range(3000):
        m = rng.choice([8, 12, 30, 73, 500, 10 ** 4])
        r0 = rng.choice([1, 1, 2, 3])
        if r0 >= m:
            continue
        nu = mpf(rng.random()) * rng.choice([0.3, 1, 3])
        mu = mpf(rng.random()) * 3 + mpf('0.01')
        rest = m - nu - mu
        if rest <= 0:
            continue
        hist = []
        # split mu over ranks 1..r0 and rest over random ranks in [1, m]
        parts = [rng.random() for _ in range(r0)]
        for r, p in zip(range(1, r0 + 1), parts):
            hist.append((r, mu * p / sum(parts) / r))
        parts = [rng.random() for _ in range(3)]
        ranks = [rng.randint(1, m) for _ in range(3)]
        for r, p in zip(ranks, parts):
            hist.append((r, rest * p / sum(parts) / r))
        a = root(hist, m)
        b = log(1 + nu / mu) / log(mpf(m) / r0)
        k += 1
        if a > b + mpf(10) ** -40:
            print('P-prime VIOLATION', m, r0, nu, mu, hist, a, b); viol += 1
        worst = max(worst, a / b)
    print(f"P' random placement check: {k} histograms, max root/bound = {float(worst):.4f}")
    for mu in [1, 2]:
        for m in [73, 10 ** 6]:
            print(f'  nu=5/14 u=1 mu={mu} m={m}: sharpened bound {math.log(1 + (5 / 14) / mu) / math.log(m):.4f}')
    print('PASS_EXPQQ_LEMMAP' if viol == 0 and n == 284 and k > 0 else 'FAIL_EXPQQ_LEMMAP')


if __name__ == '__main__':
    main()
