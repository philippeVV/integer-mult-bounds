#!/usr/bin/env python3
"""expL6 planning arithmetic (stdlib). Reuses ../../peel-dead/hyp_drop.py (expQH) functions unchanged.
Run: python3 -I -B estimates.py   (a few seconds).  Prints PASS_EXPL6_ESTIMATES."""
import math, os, sys, importlib.util
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("hyp_drop", os.path.join(here, "..", "..", "peel-dead", "hyp_drop.py"))
hd = importlib.util.module_from_spec(spec); spec.loader.exec_module(hd)
K0 = 0.0017
kap = hd.kap

def sec(t): print("\n== " + t)

def kd(delta, L):  # coding-seed conversion used by expNC/NI/QI
    y = -math.log(1 - delta); return y / (L - y)

def a_qk(k, m, b):  # expQC/QK: a = 1 - ln(m/k)/ln b
    return 1 - math.log(m / k) / math.log(b)

sec("O1/O2: outer conversion of a supplier saving a (expQH part E + assembly-domain-audit)")
print(" a          unfused a/(1+a)  +Gauss r<1/4   +guard 1/33   fused, r_max=1/2   fused, no range cap")
rows = [("today", 0.0016755), ("frame. realistic", 0.03), ("frame. formal", 0.1230),
        ("Tanner q=64 (QK)", 0.197), ("Tanner q->inf", 0.317), ("peel nu=5/14", 0.357), ("a=1", 1.0)]
for lbl, a in rows:
    u = kap(a)
    print(f" {lbl:17s} a={a:.4f}  {u:.4f} (x{u/K0:6.1f})  {min(u,.25):.4f} (x{min(u,.25)/K0:5.1f})"
          f"  {min(u,1/33):.4f}   {min(a,.5):.4f} (x{min(a,.5)/K0:5.1f})   {a:.4f} (x{a/K0:5.1f})")
print(" gain of fusing alone = (1+a): today x%.4f, at a=0.03 x%.3f, at a=0.317 x%.3f" % (1.0016755, 1.03, 1.317))

sec("I1: endpoint twists absorbed by the consumer -> drop I6 (u<=h/2) [hyp_drop sup()]")
for lbl, kw in [("all kept", {}), ("u<=h-1", dict(umax=lambda h: h - 1)), ("u<=h", dict(umax=lambda h: h)),
                ("u<=h, gamma=1/2", dict(gamma=0.5, umax=lambda h: h)),
                ("u<=20h (labels beyond stage space)", dict(umax=lambda h: 20 * h, hmax=40)),
                ("u<=20h, gamma=1/2", dict(gamma=0.5, umax=lambda h: 20 * h, hmax=40))]:
    a0, s, h, u, m = hd.sup(**kw); a = min(a0, 1.0)
    print(f" {lbl:36s} s={s:2d} h={h:3d} u={u:4d} m={m:5d} a<={a:.4f} kappa<={kap(a):.4f} (x{kap(a)/K0:6.1f}) fused x{a/K0:6.1f}")
# one-block both-bank envelope (expQT W2) at P112: gamma_p=gamma_q=64/112, h=9
g = 64 / 112
for H in (0.0, 1.0):
    d = 2 * (2 - 2 * g - g) / ((2 + H) * 9); y = -math.log(1 - d)
    print(f" QT W2 envelope P112 H/V={H}: delta<={d:.4f} kappa<=y/(2.93+y)={y/(2.93+y):.4f} (x{y/(2.93+y)/K0:.1f})")
# break-even sink-freedom: with twisted sinks the window cost could drop from gamma to gamma' ; nu = u(1 - s/(s-1) gamma')
for gp in (4/7, 0.5, 0.4, 0.3, 0.2, 0.0):
    a0, s, h, u, m = hd.sup(gamma=gp); a = min(a0, 1)
    print(f"   line-port template, gamma={gp:.3f}: a<={a:.4f} kappa<={kap(a):.4f} (x{kap(a)/K0:.1f})")

sec("C1/F1: coding-seed conversions (expNC/QI kd, expQC/QK a formula)")
for lbl, d in [("PG(2,16)-M flow cap", 0.0684), ("schedule784 exact (QE)", 0.0857), ("PG(2,16) odd (non-circuit)", 0.1687)]:
    print(f" {lbl:28s} delta={d:.4f}: kappa L=ln5 {kd(d, math.log(5)):.4f} (x{kd(d, math.log(5))/K0:.1f}),"
          f" L=ln8 {kd(d, math.log(8)):.4f} (x{kd(d, math.log(8))/K0:.1f}), L=4 {kd(d, 4):.4f} (x{kd(d, 4)/K0:.1f})")
for q, S, E in [(16, 772, 4368), (64, 57360, 266304)]:
    a = a_qk(S, E - S, 5)
    print(f" Hermitian Singer q={q}: S={S} E={E} a={a:.4f} kappa={kap(a):.4f} (x{kap(a)/K0:.1f}); fused x{a/K0:.1f}")
# char-0 (complex) twin of the q=16 layout: NK Teichmueller-lift upper bound 141 sessions
a = a_qk(141, 4368 - 141, 5)
print(f" char-0 twin of q=16 layout, S<=141 (NK): a<={a:.4f}  (negative => no complex twin from this layout)")
print(f" break-even for any square-gate girth-6 seed: a>0 iff S/E > 1/(b+1) = {1/6:.4f}; QK limit S/E -> 1/4")
for SE in (0.2, 0.25, 0.3, 1/3, 0.4, 0.45):
    print(f"   S/E={SE:.3f}, b=5: a={a_qk(SE, 1 - SE, 5):.4f}; b=11 (hexagon): a={a_qk(SE, 1 - SE, 11):.4f}")
# GQ interim (expQV q16_a.log): dim W 11332..11450 of E=73728
for S in (11332, 11450):
    a = a_qk(S, 73728 - S, 7)
    print(f" GQ T2*(O) q=16/F256 interim S={S}: S/E={S/73728:.4f}, b=7: a={a:.4f} (QV log reports 0.123-0.130)")

sec("kappa = min over suppliers: a bit-only gain")
kb, kc = 0.00173343, 0.00167546
print(f" bit-only gain, complex unchanged: kappa <= {kc:.6f} (x{kc/0.0016727:.4f} over selected, x{kc/K0:.4f} over 0.0017)")
print(f" with complex at its framework realistic cap a=0.03: kappa <= {kap(0.03):.4f} (x{kap(0.03)/K0:.1f})")

sec("kappa -> 1 bookkeeping")
print(" needs: fused (O1) + no range cap (O1/O2) + a_bit -> 1 AND a_complex -> 1.")
print(" Lemma P' (expQQ): a <= ln(1+nu/mu)/ln(m/r0). a>=0.9 at m=73,r0=1 needs nu/mu >=",
      round(73 ** 0.9 - 1, 1))
print("PASS_EXPL6_ESTIMATES")
