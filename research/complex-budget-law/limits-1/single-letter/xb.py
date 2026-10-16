"""expXB: level-scaling seeds. Stdlib only.

Part A  varying-seed supplier recursion on a dyadic log grid:
        T[i] = sum_k w_k(i) T[i-k] + h(i),  T[i<=0] = 1,  scale x = 2^i.
        Checks ln T[i] ~ sum_{j<=i} p(j) ln2  (log-scale Cesaro law), p = 1 - a.
Part B  additive per-node overhead h(i) = 2^{theta i}: harmless iff theta < p.
Part C  multiplicative per-level overhead K: a_j = (ln(1/rho) - ln K)/ln b.
Part D  amplifier comparison (repo recursion vs BGS / ZLL) on the BH seed.
Part E  ceilings: Hermitian Singer a(q); Tanner a(s, g); required s for a target;
        char-0 flat-leverage cap (expRF Markov atoms) vs requirement.
"""
import math

LN2 = math.log(2)


def root_p(w):
    """p with sum_k w[k] 2^{-k p} = 1 (w: dict shift->weight). Bisection."""
    f = lambda p: sum(c * 2.0 ** (-k * p) for k, c in w.items()) - 1
    lo, hi = 0.0, 1.0
    assert f(lo) > 0 > f(hi), (f(lo), f(hi))
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    return (lo + hi) / 2


def seed_for(a, b_shift):
    """single-rank seed: shrink 2^b_shift, weight chosen so root a (p = 1-a)."""
    p = 1 - a
    return {b_shift: 2.0 ** (b_shift * p)}


def two_rank_seed(a, k1=1, k2=3, frac=0.5):
    """two-rank seed with root exactly a: weights c1, c2 with c1 2^{-k1 p} = frac."""
    p = 1 - a
    return {k1: frac * 2.0 ** (k1 * p), k2: (1 - frac) * 2.0 ** (k2 * p)}


def run(seed_at, h, N):
    """h(i) = ln of additive overhead. log T on grid, computed in log-space to avoid overflow. returns lnT list."""
    lnT = [0.0] * (N + 1)
    for i in range(1, N + 1):
        w = seed_at(i)
        terms = [math.log(c) + (lnT[i - k] if i - k > 0 else 0.0) for k, c in w.items()]
        terms.append(h(i))  # h returns ln(overhead)
        mx = max(terms)
        lnT[i] = mx + math.log(sum(math.exp(t - mx) for t in terms))
    return lnT


def partA():
    print("== A: log-scale Cesaro law for varying seeds")
    N = 4000
    cases = {
        # a(i) -> 0.3 like 0.3 - 0.2/sqrt(i)  (Hermitian-like q^{-1/2} convergence)
        "a=0.3-0.2/sqrt(i), 2-rank": lambda i: two_rank_seed(max(0.01, 0.3 - 0.2 / math.sqrt(i))),
        # alternating seeds 0.05 / 0.4, constant shift: Cesaro mean = 0.225
        "alternate a=0.05/0.40 per 50 steps": lambda i: two_rank_seed(0.05 if (i // 50) % 2 else 0.40),
        # growing per-level shrink (b grows): a = 1 - ln(3)/ln(b), b = 2^k, k grows ~ sqrt(i)
        "Tanner-like a=1-ln3/ln b, b=2^{floor(sqrt i)+2}": lambda i: seed_for(
            1 - math.log(3) / ((int(math.sqrt(i)) + 2) * LN2), int(math.sqrt(i)) + 2),
    }
    for name, seed_at in cases.items():
        lnT = run(seed_at, lambda i: 0.0, N)
        for n in (500, 1000, 2000, 4000):
            pred = sum(root_p(seed_at(j)) for j in range(1, n + 1)) * LN2
            ces_a = 1 - pred / (n * LN2)
            eff_a = 1 - lnT[n] / (n * LN2)
            print(f"  {name:48s} i={n:5d}  a_eff(true)={eff_a:.4f}  a_Cesaro={ces_a:.4f}  "
                  f"lnT-pred={lnT[n]-pred:+.3f}")


def partB():
    print("== B: additive overhead h(i)=2^{theta i} at fixed a=0.3 (p=0.7)")
    N = 3000
    for theta in (0.0, 0.3, 0.6, 0.69, 0.72, 0.9):
        lnT = run(lambda i: two_rank_seed(0.3), lambda i, t=theta: t * i * LN2, N)
        eff = [1 - lnT[n] / (n * LN2) for n in (1000, 2000, 3000)]
        print(f"  theta={theta:4.2f}: a_eff at i=1000/2000/3000 = " + " / ".join(f"{e:.4f}" for e in eff)
              + f"   (prediction: {'0.3' if theta < 0.7 else f'{1-theta:.2f}'})")
    print("  joint growth b_j ~ exp(sqrt(C ln y)) with h ~ ln N ~ b: see README §1.4 (analytic)")


def partC():
    print("== C: multiplicative per-level overhead K (Lemma C dilution, decoder etc.)")
    rho = 147 / 157
    for b in (5, 11, 101, 10 ** 4, 10 ** 8):
        for K in (1, 3, 5):
            a = (math.log(1 / rho) - math.log(K)) / math.log(b)
            print(f"  BH-like rho, b={b:>9}: K={K}  a={a:+.4f}", end="")
        print()
    print("  Tanner s=1/4, any constant K: a = 1 - ln(3K)/ln b -> 1 as b -> inf (K amortised)")


def partD():
    print("== D: amplifiers on the Braverman-He seed (rho = 147/157, b = 5, c1 = 67.94)")
    rho, b, c1 = 147 / 157, 5, 67.94
    eps = 1 / rho - 1
    a_repo = math.log(1 / rho) / math.log(b)
    k_repo = a_repo / (1 + a_repo)
    c2 = math.log(1 + eps / 2) / math.log(16 * c1 ** 2)
    c2opt = 0.00412
    allN = c2 / 2
    G = b + 1  # ZLL Gamma >= max(3 m0, b0)+1 >= b+1 (very generous lower bound)
    zll = math.log(1 / rho) / math.log(4 * G ** 2)
    print(f"  repo recursion a = ln(1/rho)/ln b      = {a_repo:.5f}   (kappa=a/(1+a) {k_repo:.5f})")
    print(f"  BGS printed c2 = ln(1+eps/2)/ln(16c1^2) = {c2:.5f}   ratio a/c2 = {a_repo/c2:.1f}")
    print(f"  BGS alpha-optimised (expRL)            = {c2opt:.5f}   ratio = {a_repo/c2opt:.1f}")
    print(f"  BGS all-N loose c2/2 (expRL App.F: 0.00198) = {allN:.5f}   ratio = {a_repo/allN:.1f}")
    print(f"  ZLL-type ln(1/rho)/ln(4 Gamma^2), Gamma=b+1 (generous) = {zll:.5f}  ratio = {a_repo/zll:.2f}")
    print("  loss decomposition repo vs BGS: numerator ln(1/rho)/ln(1+eps/2) = "
          f"{math.log(1/rho)/math.log(1+eps/2):.2f};  denominator ln(16c1^2)/ln b = "
          f"{math.log(16*c1**2)/math.log(b):.2f}")


def a_tanner(s, g):
    return 1 - math.log((1 - s) / s) / math.log(g - 1)


def s_req(a, g):
    return 1 / (1 + (g - 1) ** (1 - a))


def partE():
    print("== E1: Hermitian Singer PG(2,q)-M (char 2, b=5): S = q(q^2-q^1.5+1)/4, E = q(q^2+q+1)")
    for q in (16, 64, 256, 1024, 4096, 2 ** 16, 2 ** 20):
        S = q * (q * q - q ** 1.5 + 1) / 4
        E = q * (q * q + q + 1)
        s = S / E
        a = a_tanner(s, 6)
        print(f"  q={q:>8}: S/E={s:.4f}  a={a:.4f}  kappa=a/(1+a)={a/(1+a):.4f}  x{a/(1+a)/0.0017:6.1f}"
              f"   a_inf-a={1-math.log(3)/math.log(5)-a:.4f}")
    ainf = 1 - math.log(3) / math.log(5)
    print(f"  a_inf = 1 - ln3/ln5 = {ainf:.4f}; kappa_inf = {ainf/(1+ainf):.4f} (x{ainf/(1+ainf)/0.0017:.0f})")
    print("== E2: Tanner a(s,g) and requirement s_req(a,g) = 1/(1+(g-1)^{1-a})")
    print("   g   a(s=1/4)  s_req(a=1/9)  s_req(a=1/2)  flat-char0 cap M_g  a_cap(flat)")
    Mg = {4: 0.25, 6: 2 / 11, 8: 0.125, 10: 0.1034, 12: 1 / 12, 16: 1 / 16}
    for g in (4, 6, 8, 10, 12, 16, 20, 32, 64):
        m = Mg.get(g)
        mtxt = f"{m:.4f}" if m else "  -   "
        acap = f"{a_tanner(m, g):+.3f}" if m else "  -  "
        print(f"  {g:3d}   {a_tanner(0.25, g):+.3f}     {s_req(1/9, g):.4f}       {s_req(0.5, g):.4f}"
              f"       {mtxt}           {acap}")
    print("  (even g/2: M_g = 1/g = s_req(a=0); kappa>=0.1 needs a>=1/9 > every flat-leverage cap)")


if __name__ == "__main__":
    partA()
    partB()
    partC()
    partD()
    partE()
    print("PASS_EXPXB")
