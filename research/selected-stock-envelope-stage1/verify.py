#!/usr/bin/env python3
"""Stage-1 selected stock envelope: source bindings + exact-integer check.

Usage: python3 -I -B research/selected-stock-envelope-stage1/verify.py [repo-root]
Exact integers / Fractions only. Reads sources; executes none of them.
"""
import hashlib, json, sys
from fractions import Fraction as Q
from pathlib import Path

R = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[2]
def need(c, m):
    if not c: raise SystemExit('FAIL ' + m)
def line(path, n, text):
    ls = (R / path).read_text().splitlines()
    need(text in ls[n - 1], f'{path}:{n} must contain {text!r}')
    return f'{path}:{n}'

# 1. Bindings: which producer is bootstrapped, and its parameters.
cites = [
 line('research/bootstrap-depth-audit/sources/pr387/refine.py', 66, 'degraded_bit_saving=str(saving),bit_coarse=str(c)'),
 line('research/split-centres-global/proofs/PR386-BIT-GEOMETRY.txt', 147, 'one-level halving is false. Two levels suffice'),
 line('research/split-centres-global/proofs/PR386-BIT-GEOMETRY.txt', 148, '2*42^2=3528 < 64^2=4096'),
 line('research/split-centres-global/proofs/PR386-BIT-GEOMETRY.txt', 163, 'bootstrap cutoff inequalities'),
 line('research/split-centres-global/proofs/PR386-BIT-GEOMETRY.txt', 166, 'inherited coarse q^(m^2+1) bound'),
 line('research/five-stage-source527-banks/proof/three-stage-cover-bit.tex', 198, 'Take $q>2^{80}$'),
 line('research/five-stage-source527-banks/proof/three-stage-cover-bit.tex', 211, 'Let $W_0=|\\operatorname{GL}_m(\\mathbb F_q)|'),
 line('research/five-stage-source527-banks/proof/three-stage-cover-bit.tex', 212, 'The full role count is $W(w)=W_0q^{m^2(w-1)}$'),
 line('research/five-stage-source527-banks/proof/three-stage-cover-bit.tex', 256, 'T(n,w)\\le\\sum_rc_rT(r\\lfloor n/m\\rfloor,w)'),
 line('research/five-stage-source527-banks/proof/three-stage-cover-rows.tex', 5, 'A sufficient common divisor is $W(w)^{D(n)}$ times the inherited old-leaf'),
 line('research/five-stage-source527-banks/proof/three-stage-cover-rows.tex', 7, 'logarithmic reserve'),
 line('research/coordinated-frames-and-entrance-banks/references/pr141/ORIGINAL-PROOF.md', 9, 'old ordinary saving 384599/10^10'),
 line('research/coordinated-frames-and-entrance-banks/references/pr141/ORIGINAL-PROOF.md', 25, 'W_complex^D_complex * W_coarse^D_coarse * W_old^D_old'),
 line('research/coordinated-frames-and-entrance-banks/references/pr141/ORIGINAL-PROOF.md', 29, 'plus a fixed constant'),
 line('research/coordinated-frames-and-entrance-banks/references/pr141/ORIGINAL-PROOF.md', 35, '| Old ordinary leaf | 575 | 529 | 28 | 9 | 252 |'),
 line('research/bootstrap-depth-audit/sources/pr387/PROOF.md', 7, 'performs required high-digit borrowing internally'),
 line('upstream/build/sections/04-swap.tex', 354, 'by each required $W^j$'),
]
cites = [c for c in cites if c]

asm_p = R / 'research/community-frontier1094-audit/ASSEMBLY.json'
sel = json.loads((R / 'certificates/selected-result.json').read_text())
need(hashlib.sha256(asm_p.read_bytes()).hexdigest() == sel['certificate_sha256'], 'selected record binds ASSEMBLY')
asm = json.loads(asm_p.read_text())
exp_p = R / 'research/e8-bit-p111-tuned/tune/EXPECTED.json'
need(asm['input_sha256']['package:inputs/PR589-EXPECTED.json'] == hashlib.sha256(exp_p.read_bytes()).hexdigest(), 'ASSEMBLY pins PR589')
bit = json.loads(exp_p.read_text())['nine']['bit']
m, rmax, depth, phys = 64, bit['maxchild'], bit['halving_depth'], bit['physical']
need((rmax, depth, phys) == (42, 2, 28425), 'selected bit producer parameters')
need(bit['row_coefficient'] == depth * (m * m + 1) + 9909 + 252 == 18355, 'external row formula binds m=64, depth 2')
need(m ** 2 > 2 * rmax ** 2 and not m > 2 * rmax, 'halving depth 2 is minimal')
need(bit['two_level_coefficient_bits'] < bit['charged_two_level_coefficient_bits'] == 512, 'cutoff 512 = bit charged coefficient')

# 2. Depth lemma D(n) <= 2*ceil(log2 max(2,n)) for the recurrence n -> r*floor(n/64), r<=42, leaf n<64.
#    Proof: two levels map n to <= (42/64)^2 n < n/2; a level exists only while size >= 64 >= 1.
def D(n):
    d = 0
    while n >= m: n, d = rmax * (n // m), d + 1
    return d
clog = lambda n: (max(2, n) - 1).bit_length()
need(all(D(n) <= depth * clog(n) for n in range(1, 1 << 20)), 'depth lemma sample')
need(2 * rmax ** 2 < m ** 2, 'two-level contraction')

# 3. Per-level stock: log_q W(w) = log_q W_0 + m^2(w-1) < m^2 w + 1 when |G_q| <= |GL_64(F_q)| < q^4096
#    and the per-level role multiplier N < q. N <= physical * (1 + 32 m^2 * paid children) is far below 2^80.
paid = sum(json.loads(exp_p.read_text())['nine']['profile'].values())
N = phys * (1 + 32 * m * m * paid)
need(N < 2 ** 80, 'role multiplier below q')
# Chain: log_q S <= D(n)(4096w+1) + log_q R0 <= 2(1+t)*4097*w + (252*(1+t)+28c)/80,  t=log2 e>=1, w>=1.
SIGMA_NEW = depth * (m * m + 1)                     # 8194
need(SIGMA_NEW == 8194, 'parent coefficient')
# sigma_1 <= 8194 + 252/80 + 28c/80  <=  8198 + c   (c = PR141 'fixed constant' depth term, c>=0)

# 4. Stage-1 guard (reproduction of recorded numbers, not a new result).
for key in ('bootstrap', 'conservative_bootstrap'):
    g, L = Q(asm[key][0]['gap']), int(asm[key][0]['cutoff_log2'])
    x = L * g / 2; H = x.numerator // x.denominator - 4
    need(H == 6562783, 'H_1')
    # Guard: log2(8 sigma) <= H  <=>  8 sigma <= 2^H (sigma integer bound).
    need(8 * SIGMA_NEW <= 1 << H, 'internal-borrowing reading passes')
    need(8 * (SIGMA_NEW + 4) <= 1 << H, 'old-leaf reading passes at c=0')
    cmax = (1 << (H - 3)) - 8198
    need(8 * (8198 + cmax) <= 1 << H and 8 * (8198 + cmax + 1) > 1 << H, 'exact c threshold')
    # Alternative if the producer were the m=73 complex cover (halving depth 1: 73^2 > 2*48^2):
    need(73 ** 2 > 2 * 48 ** 2 and 8 * (1 * (73 * 73 + 1) + 4) <= 1 << H, 'm73 alternative')
print(json.dumps(dict(status='PASS_STAGE1_CONDITIONAL_ENVELOPE', producer='PR589 bit m=64', depth=depth,
    per_level='4096w+1', sigma_new=SIGMA_NEW, sigma_internal_reading=SIGMA_NEW,
    sigma_old_leaf_reading='8198+c', H1=6562783, c_threshold_bits=6562780,
    role_multiplier_bits=N.bit_length(), citations=cites), indent=1))
