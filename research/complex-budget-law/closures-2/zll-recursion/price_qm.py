"""expQM / price_qm.py -- ideal price of a k-level unrolled word that pays a per-level overhead ONCE.

Engines (pinned, imported unchanged via engines.py (slim copy of expM8/price.py)):
  complex: research/complex-local-word-ceiling/ceiling.py  relax() = per-call moment Psi_S(a)  (m=73, 896 fam, 9 st, R=51)
  bit    : research/bit-revision3-levers-seed03/levers.py  bill(relaxed()) = Psi_S(a)          (m=64, w=888/9, R=45)
Model (expQG sec.3 + this note sec.2):
  * a complete level S has moment Psi_S; unrolling j complete levels gives Psi_S^j (tensor invariance).
  * a k-level word that pays overhead X once = S (x) S_X^(k-1), where S_X is the level with X deleted
    (its mass handled as the engine handles it: centre copies removed and charge 0; helper/centre climbs -> idle at R).
    Its moment is Psi_S * Psi_SX^(k-1) (expQG Thm (i): Psi multiplies under unrolling, ranks multiply, widths multiply).
  * the root a_k of Psi_S Psi_SX^(k-1) = 1 is a secant-weighted average of a_S and a_SX, -> a_SX as k -> inf.
  * kappa = a/(1+a), a = min(complex, bit) (two suppliers; expM7).
  * case (c) of README sec.3: only calls that run INSIDE an outer centre child skip their centre (shared source);
    two-type branching moment, root of the Perron value.
Run: python3 -I -B price_qm.py   (stdlib + pinned engines; ~1 min)
"""
import sys, os
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
from fractions import Fraction as Q
import engines as P
Cx, Bt = P.Cx, P.Bt
REF = Q('0.0017')

# ---- per-call local profiles (selected words) -----------------------------------------------------------
CXD = {8: 224}                                    # complex data floor, one-shot (favourable: largest admissible rank)
BXD = {7: 222}                                    # bit data floor one-shot (2V = 222 registers x 7)


def psi_cx(local, centre, a):
    lo, hi = Cx.relax(local, 73, 896, 9, centre, a)
    return float(lo), float(hi)


def psi_bt(local, cst, a):
    H = Bt.relaxed(local, 64, Q(888, 9), cst, 45)
    lo, hi = Bt.bill(H, 64, Q(888, 9), a)
    return float(lo), float(hi)


LEVELS = {
    # name: (complex (local, centre), bit (local, cst))
    'S (selected level)': (({**P.CW, 8: 8}, 64), ({**P.BW, 7: P.BW[7] + 9}, 63)),
    'S_C (centre/decoder deleted)': ((dict(P.CW), 0), (dict(P.BW), 0)),
    'S_CH (centre + helpers deleted; data floor + idle)': ((dict(CXD), 0), (dict(BXD), 0)),
}


def psi(name, side, a):
    (cl, cc), (bl, bc) = LEVELS[name]
    return psi_cx(cl, cc, a) if side == 'cx' else psi_bt(bl, bc, a)


def root_word(parts, side, lo=1e-6, hi=0.06, it=60):
    """root of prod_j Psi_{name_j}(a)^{e_j} = 1 (upper ends of the intervals -> the ceiling side)."""
    def f(a):
        import math
        v = 0.0
        for name, e in parts:
            v += e * math.log(psi(name, side, Q(a).limit_denominator(10 ** 15))[0])
        return math.exp(v) if v < 50 else float('inf')
    for _ in range(it):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(mid) < 1 else (lo, mid)
    return hi


def root_two_type(side, it=55):
    """case (c): calls inside an outer centre child are centre-free (type B); everything else type A.
    Perron root of [[PsiA - psi_c, psi_c], [PsiB, 0]] = 1  <=>  (PsiA - psi_c) + psi_c * PsiB = 1."""
    lo, hi = 1e-6, 0.06
    for _ in range(it):
        a = (lo + hi) / 2; qa = Q(a).limit_denominator(10 ** 15)
        A = psi('S (selected level)', side, qa)[0]
        B = psi('S_C (centre/decoder deleted)', side, qa)[0]
        if side == 'cx':
            pc = 8 * float(Cx.child_cost(8, 73, qa)[0]) / (896 / 9)
        else:
            pc = 9 * float(Bt.child(7, 64, qa)[0]) / (888 / 9)
        v = A - pc + pc * B
        lo, hi = (a, hi) if v < 1 else (lo, a)
    return hi


def kap(a):
    return a / (1 + a)


def main():
    base = {s: root_word([('S (selected level)', 1)], s) for s in ('cx', 'bt')}
    print(f"single complete level: complex a={base['cx']:.7f}  bit a={base['bt']:.7f}  kappa={kap(min(base.values())):.6f}")
    # tensor invariance control
    t2 = {s: root_word([('S (selected level)', 2)], s) for s in ('cx', 'bt')}
    assert all(abs(t2[s] - base[s]) < 1e-10 for s in t2), t2
    print('control: S(x)S root = S root on both sides (tensor invariance)')
    k_ref = kap(min(base.values()))
    rows = {}
    for X in ('S_C (centre/decoder deleted)', 'S_CH (centre + helpers deleted; data floor + idle)'):
        print(f'\nk-level word paying the overhead once; deleted at k-1 levels: {X}')
        for k in (2, 3, 4, 9, 10 ** 6):
            parts = [('S (selected level)', 1), (X, k - 1)]
            ac, ab = root_word(parts, 'cx'), root_word(parts, 'bt')
            kk = kap(min(ac, ab)); kc = kap(min(ac, base['bt']))
            rows[(X, k)] = (ac, ab)
            print(f'  k={k:>7}: complex a={ac:.6f} (x{ac / base["cx"]:.3f})  bit a={ab:.6f} (x{ab / base["bt"]:.3f})  '
                  f'kappa={kk:.6f} x{kk / float(REF):.2f} over 0.0017;  complex-only (bit fixed) x{kc / float(REF):.3f}')
    # reconciliation with the lever rows
    c_free = Cx.ceiling(dict(P.CW), 73, 896, 9, 0, grid=Q(1, 10 ** 12))[1]
    sh2 = Cx.ceiling({**P.CW, 8: 4}, 73, 896, 9, 32, grid=Q(1, 10 ** 12))[1]
    c0 = Cx.ceiling({**P.CW, 8: 8}, 73, 896, 9, 64, grid=Q(1, 10 ** 12))[1]
    print(f'\nreconcile: levers.py "centre free" x{float(c_free / c0):.3f} = k->inf limit x{rows[("S_C (centre/decoder deleted)", 10**6)][0] / float(c0):.3f};'
          f' levers.py "shared by 2 calls" (C=32, one level) x{float(sh2 / c0):.3f} vs 2-level mediant '
          f'x{rows[("S_C (centre/decoder deleted)", 2)][0] / float(c0):.3f}')
    # universal Lambda cap (budget law sec.6): u/((m-u) ln(m/R)), the x8.2 * x2.8 ~ x23 relaxation
    import math
    uni = 1 / (72 * math.log(73 / 51)); jens = math.log(73 / (73 - 5 / 14)) / math.log(73 / 51)
    print(f'universal caps (complex, all paid mass at rank R=51): u=1 a<={uni:.5f} (x{uni / float(c0):.1f}); '
          f'at nu=5/14 (Lambda-only lever) a<={jens:.6f} (x{jens / float(c0):.2f})')
    acx, abt = root_two_type('cx'), root_two_type('bt')
    kk = kap(min(acx, abt))
    print(f'\ncase (c): only calls inside outer centre children skip their centre (shared source), all depths: '
          f'complex x{acx / base["cx"]:.4f}  bit x{abt / base["bt"]:.4f}  kappa x{kk / k_ref:.4f} over selected')
    print('PASS_EXPQM_PRICE')


if __name__ == '__main__':
    main()
