"""E1 numerics: gap exponents guaranteed by the published refutations, vs kappa.
BGS17 Thm 1.1 / App. F: from a seed with gap 1+eps the family has gap
Omega((log size)^{c2}), c2 = log(1+eps/2)/log(16 c1^2), c1 = l/w + 1 >= 2.
ZLL (2610.09367) eq. (31): c = log(1/rho0)/log(4 Gamma^2), Gamma >= b0+1 > 2*390625.
"""
import math
eps = 157/147 - 1                      # Braverman-He seed: coding 1 vs flow <= 147/157
print(f"BH seed gap 1+eps, eps = {eps:.5f}")
for c1 in (2, 3, 5, 10, 30, 100):
    c2 = math.log(1+eps/2)/math.log(16*c1*c1)
    print(f"  c1={c1:4d}: BGS exponent c2 = {c2:.5f}  c2/kappa(0.0017) = {c2/0.0017:.2f}")
lnLam_min = math.log(4*(2*390625+1)**2)
print(f"ZLL: log Lambda >= {lnLam_min:.2f}, so c <= {1/lnLam_min:.4f} * ln(1/rho0)")
# AFKL accounting: oblivious T = n lg^{1-kappa} n  =>  gamma >= lg n /(4 C log2(2c) lg^{1-kappa} n)
for k in (0.0017, 0.03, 0.317, 1.0):
    for lgn in (64, 2**20):
        print(f"  kappa={k}: lg n={lgn}: lg^kappa n = {lgn**k:.4g}")
