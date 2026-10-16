#!/usr/bin/env python3
"""expQT pricing for idea 8 (endpoint labels on both banks). Stdlib, < 1 s.
Run: python3 -I -B price.py  -> prints PASS_EXPQT_PRICE.
Every number is a favourable planning value (upper envelope), not a claim."""
import math

BASE = 0.0017
LBAR = 2.93
def kap(delta, L=LBAR):
    if delta <= 0: return 0.0
    y = -math.log(1 - delta); return y / (L + y)

print("== §1 one block (Theorem W vs Theorem W2), h = 9, W = 2V + H helpers, kappa by gadget law at lbar = 2.93")
h = 9
rows = []
for name, gp, gq in (("P120 both E8", 72/120, 72/120), ("P112 both E8", 64/112, 64/112),
                     ("gamma 1/2", 0.5, 0.5), ("gamma 0.4", 0.4, 0.4)):
    one = (1 - 2 * gp) / h                      # Thm W (one bank, data wires never drop)
    for k in (0, 1, 8):                         # helper wires per port (k = 8: one per Gram term)
        both = 2 * (2 - 2 * gp - gq) / (h * (2 + k))   # Thm W2: D = 2(2V - 2C_p - C_q), Wm = (2+k)V h
        rows.append((name, k, one, both))
        print(f"  {name:13s} H={k}V: one-bank delta<={one:+.4f}   both-bank delta<={both:+.4f}"
              f"  kappa<={kap(both):.4f} (x{kap(both)/BASE:5.2f})")
# threshold: both-bank one-block column open iff 2*gp + gq < 2  (gp = gq: gamma < 2/3)
assert abs(rows[3][3] - 2 * (2 - 3 * 64/112) / (9 * 2)) < 1e-12

print("== §2 template (expL5 §D formula): nu' = nu + u(1 - tau' g), s = 9, h = 9, u = 1, g = 4/7")
s, u, g = 9, 1, 4/7; tau = s / (s - 1); m = u + s * (h - 1)
nu = u * (1 - tau * g); e0 = nu / m
be = (1 - (nu * (m + u) / m - nu) / u) / g
print(f"  today eps = {e0:.5f}; break-even tau' (m' = m + u) = {be:.3f}")
for lab, mp in (("m' = m + u  (expL5)", m + u), ("m' = m + h - 1 (one more call)", m + h - 1)):
    for taup in (1.0, tau, 1.44, be):
        nup = nu + u * (1 - taup * g); e1 = nup / mp
        print(f"  {lab:32s} tau'={taup:.3f}: eps {e1:.5f}  x{e1/e0:.2f}")
print(f"  Thm W2 forces only tau' >= tau * gamma_q / gamma_p = {tau:.3f} for a congruent Y system (gamma_q = gamma_p)")

print("== §3 coverage law (expSB / #1587 follow-up 5): per-call deficit = (8/9) tau_unpaid - C_total")
inc = (8 / 9) * 112 - 64
for lab, tu, C in (("incumbent lines", 112, 64), ("both banks, two Gram decoders", 224, 128),
                   ("both banks, Y window private (C_q = 112)", 224, 64 + 112)):
    dd = (8 / 9) * tu - C
    print(f"  {lab:42s} tau={tu} C={C}: deficit/call {dd:7.2f}  (x{dd/inc:.2f} of incumbent)")
print("  BUT the rank invariant (#1587 c18914510 s2, U <= v(N-d), #1094 at equality) forbids tau_unpaid = 224 at N = 9")
print("  with line labels; the extra unpaid coordinates need more calls: ctoth's priced two-unpaid template = -15.5% (x0.845)")

print("== §4 suppliers: kappa = min(complex, bit); bit endpoint caps complex-only gains at +1.1%")
for lab, xc, xb in (("template (invariant)", 0.845, 1.0), ("template, expL5 optimistic (infeasible)", 1.97, 1.0),
                    ("one-block W2 envelope P112, H=V, both suppliers", kap(rows[4][3]) / BASE, None)):
    tot = min(xc, 1.011) if xb == 1.0 else xc
    print(f"  {lab:48s} complex x{xc:.2f}  bit {'x1.00 (oracle: Y labels kill the motif deficit)' if xb else 'needs a bit word'}  -> overall x{tot:.3f}")
print("PASS_EXPQT_PRICE")
