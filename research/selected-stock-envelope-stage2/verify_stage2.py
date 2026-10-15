#!/usr/bin/env python3
"""Selected stock envelope at native bootstrap stage 2 (card selected-stock-envelope, rev. 4).

Usage: python3 -I -B verify_stage2.py <repo-root>
Exact integers / Fractions only. Reads repository files; executes none of them.

Stage 2 is the first native stage whose stopped parent calls a *stopped* completed leaf
(A_1) rather than the base leaf A_0. This checker binds the stage-2 parameters to the
selected record, inventories the source-stated widths of the A_1 calls made by the
stage-2 body, and checks the exact inequalities of the stage-2 envelope

    log_q S_2(e,w) <= (8194 + eps_2(gamma)) * w * (1 + log2 e)    for all e >= 2^L_2,

under each of the three readings of the completed inner leaf (internal / retained /
fully hoisted), with gamma = log2 of the routing-width constant c_R. See README.md.
"""
import hashlib, json, sys
from fractions import Fraction as Q
from pathlib import Path

R = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path.cwd()


class Fail(Exception):
    pass


def need(c, m):
    if not c:
        raise Fail(m)


def line(path, n, text):
    ls = (R / path).read_bytes().decode().split('\n')
    need(n <= len(ls) and text in ls[n - 1], f'{path}:{n} must contain {text!r}')
    return f'{path}:{n}'


PR185 = 'research/community-round8-audit/pr185-proof-reviewed.md'
PR387 = 'research/bootstrap-depth-audit/sources/pr387/PROOF.md'
ROWS = 'research/five-stage-source527-banks/proof/three-stage-cover-rows.tex'
BIT = 'research/five-stage-source527-banks/proof/three-stage-cover-bit.tex'
PR141 = 'research/coordinated-frames-and-entrance-banks/references/pr141/ORIGINAL-PROOF.md'
CITES = [
    # Which ordinary calls the stage-j stopped parent makes, and at which widths.
    (PR387, 9, 'using only completed A_(j-1) for ordinary leaves, routing, peeled atoms, bounded stopped leaves and high-part exchanges'),
    (PR185, 13, 'Weighted leaves: four ordinary atom interchanges'),
    (PR185, 14, 'Peeled atoms and bounded stopped leaves: a fixed number of ordinary width-`w` calls'),
    (PR185, 15, 'Fixed cover routing: a fixed number of ordinary width-`O(w)` calls'),
    (BIT, 208, 'a constant number of old width-$O(w)$ interchanges'),
    (BIT, 258, 'T(n,w)\\le Cw^{1-a_{\\rm old}}\\quad(n<m)'),
    (BIT, 260, 'Peeling fewer than $m$ leftover atoms uses the weighted leaf algorithm.'),
    (PR185, 40, 'The parent high-part exchange has width `k=O(w log e)`'),
    (ROWS, 44, 'Execute the completed \\emph{ordinary} body wrapper with this row reserve,'),
    (ROWS, 50, 'One inherited ordinary width-$k$ interchange of the two high parts,'),
    (ROWS, 65, 'its complete atom body and the old ordinary supplier to fewer than $w$'),
    # Reserve composition: one inherited leaf reserve, reused; inner borrowing readings.
    (ROWS, 5, 'A sufficient common divisor is $W(w)^{D(n)}$ times the inherited old-leaf'),
    (ROWS, 14, 'so successive residuals and the outer wrapper reuse one reserve instead'),
    (PR185, 38, 'its new ordinary leaves need no *external* extra stock'),
    (PR185, 38, 'it borrows from the two chunks passed to that ordinary call'),
    (PR387, 7, 'performs required high-digit borrowing internally'),
    (PR387, 7, 'It is retained conditionally, not proved anew here.'),
    # Parent per-level divisor and the A_0 row.
    (BIT, 198, 'Take $q>2^{80}$'),
    (BIT, 212, 'The full role count is $W(w)=W_0q^{m^2(w-1)}$'),
    (PR141, 35, '| Old ordinary leaf | 575 | 529 | 28 | 9 | 252 |'),
    ('research/bootstrap-depth-audit/sources/pr387/refine.py', 66, 'degraded_bit_saving=str(saving),bit_coarse=str(c)'),
]

A0 = Q(384599, 10**10)
BETA = Q(1, 10**24)
clog = lambda n: (max(2, n) - 1).bit_length()          # ceil(log2 max(2,n))


def chain(bit, cx):
    s = min(bit, (1 - BETA) * cx) - Q(1, 10**27)
    a, rows = A0, []
    for _ in range(7):
        old, a = a, (1 - s) * s + s * a
        rows.append((a, min(s - a, 1 - a - s, 1 - a - s * (1 - old), 1 - s)))
    return s, rows


def ceilq(x):
    return -((-x.numerator) // x.denominator)


def floorq(x):
    return x.numerator // x.denominator


def depth(m, r, n):
    d = 0
    while n >= m:
        n, d = r * (n // m), d + 1
    return d


def stage_params(asm, key, side):
    s, rows = chain(Q(asm['bit_lower']), Q(asm[side]))
    out = []
    for j, ((a, g), rec) in enumerate(zip(rows, asm[key]), 1):
        need(Q(rec['exponent']) == a and Q(rec['gap']) == g, f'{key}[{j}] exponent/gap from a0')
        L = int(rec['cutoff_log2'])
        need(L == max(ceilq(36 / g**2), ceilq(Q(2 * (4 + 512)) / g)), f'{key}[{j}] recorded cutoff formula')
        H = floorq(L * g / 2) - 4
        need(L * g * g >= 36 and L * g >= 2 * (4 + H) and not L * g >= 2 * (4 + H + 1), f'{key}[{j}] headroom')
        out.append(dict(a=a, delta=g, L=L, H=H))
    return s, out


# ---- the stage-2 envelope -------------------------------------------------------------
LOG2_RHO1_CEIL = 14          # rho_1 = 8194 + 63/10 < 2^14, checked below


def gamma_max(s, L1, L2):
    """Largest integer gamma with (a) gamma+1 <= (1-s)L_1 and (b) eps_2(gamma) <= 1, where
    log2 eps_2(gamma) <= LOG2_RHO1_CEIL + 1 + s(gamma+1) - s(1-s)L_2."""
    ga = floorq((1 - s) * L1) - 1
    gb = floorq((s * (1 - s) * L2 - LOG2_RHO1_CEIL - 1) / s) - 1
    return min(ga, gb), ga, gb


def eps_log2_upper(s, L2, gamma):
    """Rational upper bound on log2 eps_2(gamma)."""
    return LOG2_RHO1_CEIL + 1 + s * (gamma + 1) - s * (1 - s) * L2


def envelope_ok(s, L1, L2, gamma, sigma_par, rho1):
    """All side conditions used in the written proof, as exact rational inequalities."""
    return (gamma >= 0
            and rho1 < 2 ** LOG2_RHO1_CEIL
            and gamma + 1 <= (1 - s) * L1          # A_1's own body calls fit inside its width
            and gamma + 1 <= (1 - s) * L2          # 1+log2(c_R w) <= 1+log2 e
            and s * L2 >= L1                       # width-w calls of stage 2 hit A_1's stopped branch
            and eps_log2_upper(s, L2, gamma) <= 0)  # eps_2 <= 1


def guard(sigma, H):
    """Sufficient: log2(8 sigma) <= H (as in PR1965/PR1990)."""
    return ceilq(8 * Q(sigma)).bit_length() <= H


def main():
    cites = [line(*c) for c in CITES]
    asm_p = R / 'research/community-frontier1094-audit/ASSEMBLY.json'
    sel = json.loads((R / 'certificates/selected-result.json').read_text())
    need(hashlib.sha256(asm_p.read_bytes()).hexdigest() == sel['certificate_sha256'], 'selected record binds ASSEMBLY')
    asm = json.loads(asm_p.read_text())
    need('M=64' in asm['bit_supplier']['deficit_formula'] and asm['bit_supplier']['pr'] == 589, 'ASSEMBLY bit supplier PR589, M=64')
    exp_p = R / 'research/e8-bit-p111-tuned/tune/EXPECTED.json'
    need(asm['input_sha256']['package:inputs/PR589-EXPECTED.json'] == hashlib.sha256(exp_p.read_bytes()).hexdigest(), 'ASSEMBLY pins PR589 record')
    bit = json.loads(exp_p.read_text())['nine']['bit']
    need((bit['live'], bit['physical'], bit['maxchild'], bit['halving_depth'], bit['row_coefficient'])
         == (28416, 28425, 42, 2, 18355), 'PR589 producer fields')
    cr = json.loads((R / 'research/copied-reversed/certificate.json').read_text())
    need(cr['prior']['bit'] == '384599/10000000000', 'A0 saving')
    fb = cr['finite_bridge']['bit']
    need((fb['m'], fb['maxchild'], fb['wire_bits'], fb['halving_degree']) == (575, 529, 28, 9), 'A0 m575 fields')

    # (P) parent: log_q W_q(w) < 4096 w + 1 <= 4097 w, since 28425 < 2^80 < q and |GL_64(F_q)| < q^4096;
    # D(n) <= 2 ceil(log2 max(2,n)) because 2*42^2 < 64^2 (exhaustive sanity range below).
    m, r = 64, bit['maxchild']
    need(bit['physical'] < 2 ** 80 and 2 * r * r < m * m, 'parent divisor/halving')
    need(all(depth(m, r, n) <= 2 * clog(n) for n in range(1 << 20)), 'BIT depth lemma on [0,2^20)')
    SIGMA_PAR = 2 * (m * m + 1)
    need(SIGMA_PAR == 8194 and bit['row_coefficient'] == SIGMA_PAR + 9909 + 252, 'sigma_par = 8194, matches 18355 row')
    # (R0) A_0 reserve: log2 R_0 <= 28 * 9 * ceil(log2 max(2,x)); log_q R_0 < (63/20)(1+log2 x).
    need(575 ** 9 > 2 * 529 ** 9 and not 575 ** 8 > 2 * 529 ** 8 and 28 * 9 == 252, 'A0 halving degree 9, coefficient 252')
    RHO0 = Q(252, 80)
    need(RHO0 == Q(63, 20), 'rho0')
    RHO1 = SIGMA_PAR + 2 * RHO0                 # fully hoisted A_1 external need coefficient
    need(RHO1 == Q(82003, 10) and RHO1 < 2 ** LOG2_RHO1_CEIL and RHO1 >= 2 ** (LOG2_RHO1_CEIL - 1), 'rho1 bits')

    s, nat = stage_params(asm, 'bootstrap', 'complex_lower')
    s_c, con = stage_params(asm, 'conservative_bootstrap', 'conservative_complex_lower')
    st1, st2 = nat[0], nat[1]
    need(st1['L'] == 4785575287428 and st2['L'] == 1704765459304159860, 'native L_1, L_2')
    need(st1['H'] == 6562783, 'H_1')

    # Branch lemma: at stage 1 every width-w call is below... irrelevant (A_0 is the only leaf);
    # at stage 2, w >= 2^(s L_2) >= 2^(L_1): width-w and width-k calls use A_1's stopped branch.
    need(s * st2['L'] >= st1['L'], 'stage-2 atom width exceeds stage-1 cutoff')
    branch_margin = floorq(s * st2['L']) - st1['L']

    gmax, ga, gb = gamma_max(s, st1['L'], st2['L'])
    results = {}
    for gamma in (0, 12, 18, 64, gmax):
        need(envelope_ok(s, st1['L'], st2['L'], gamma, SIGMA_PAR, RHO1), f'envelope side conditions at gamma={gamma}')
        lam = floorq(-eps_log2_upper(s, st2['L'], gamma))
        sig = SIGMA_PAR + 1                      # 8194 + eps_2, eps_2 <= 2^-lam <= 1
        need(guard(sig, st2['H']), 'stage-2 guard')
        results[str(gamma)] = dict(log2_eps2_at_most=-lam, sigma2_upper=str(sig))
    # gamma_max is sharp for the stated side conditions
    need(not envelope_ok(s, st1['L'], st2['L'], gmax + 1, SIGMA_PAR, RHO1), 'gamma_max + 1 rejected')

    # Side check (not the scoped claim): the conservative stage 2 meets the same side conditions.
    gmax_c = gamma_max(s_c, con[0]['L'], con[1]['L'])[0]
    need(all(envelope_ok(s_c, con[0]['L'], con[1]['L'], g, SIGMA_PAR, RHO1) for g in (0, 12, 64, gmax_c))
         and guard(SIGMA_PAR + 1, con[1]['H']), 'conservative stage 2 side check')

    # Fallback without the routing-width lemma (only x <= e): hoisted reading sigma_2 <= 8194 + rho_1.
    sig_fallback_H = SIGMA_PAR + RHO1
    need(sig_fallback_H == Q(163943, 10) and guard(sig_fallback_H, st2['H']), 'fallback guard')

    # ---- mutation controls: each must be rejected -----------------------------------------
    controls = {}
    def ctl(name, f):
        try:
            f()
        except Fail:
            controls[name] = 'rejected'
            return
        raise Fail(f'control {name} was accepted')
    ctl('stage1_cutoff_as_stage2', lambda: need(s * st1['L'] >= st1['L'], 'x'))
    ctl('gamma_above_max', lambda: need(envelope_ok(s, st1['L'], st2['L'], gmax + 1, SIGMA_PAR, RHO1), 'x'))
    ctl('rho1_one_bit_short', lambda: need(envelope_ok(s, st1['L'], st2['L'], 12, SIGMA_PAR, Q(2 ** 14)), 'x'))
    ctl('L2_minus_one', lambda: need(st2['L'] - 1 == max(ceilq(36 / st2['delta']**2), ceilq(Q(1032) / st2['delta'])), 'x'))
    ctl('false_halving_child46', lambda: need(2 * 46 * 46 < 64 * 64, 'x'))
    ctl('wrong_physical', lambda: need(28426 == bit['physical'], 'x'))
    ctl('A0_degree_8', lambda: need(575 ** 8 > 2 * 529 ** 8, 'x'))
    ctl('stage2_guard_at_H2_plus', lambda: need(guard(2 ** (st2['H'] - 2), st2['H']), 'x'))
    ctl('cite_wrong_line', lambda: line(PR185, 16, 'Fixed cover routing'))
    ctl('import_57604', lambda: need(SIGMA_PAR == 57604, 'x'))

    return dict(
        status='PASS_SELECTED_STOCK_STAGE2_CONDITIONAL_ENVELOPE',
        stage='native bootstrap stage 2 (ASSEMBLY.json bootstrap[1])',
        s=f'{s.numerator}/{s.denominator}', s_decimal=f'{float(s):.18f}',
        L1=st1['L'], L2=st2['L'], H2=st2['H'],
        conservative_L2=con[1]['L'], conservative_H2=con[1]['H'], conservative_gamma_max=gmax_c,
        sL2_minus_L1_floor=branch_margin,
        sigma_par=SIGMA_PAR, rho0='63/20', rho1='82003/10',
        gamma_max=gmax, gamma_max_limits=dict(A1_width=ga, eps_le_1=gb),
        envelope_by_gamma=results,
        sigma2=dict(internal='8194', retained='8194 + eps', hoisted='8194 + eps',
                    hoisted_without_width_lemma='163943/10'),
        guard_bits_needed=ceilq(8 * Q(SIGMA_PAR + 1)).bit_length(),
        citations=len(cites), controls=controls)


if __name__ == '__main__':
    try:
        res = main()
    except Fail as exc:
        print(json.dumps(dict(status='FAIL', reason=str(exc))))
        sys.exit(1)
    print(json.dumps(res, indent=1))
