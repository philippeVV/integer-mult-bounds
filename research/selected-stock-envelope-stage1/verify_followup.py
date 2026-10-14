#!/usr/bin/env python3
"""Selected stock envelope follow-up: old-leaf constant c and all 14 recorded stages.

Usage: python3 -I -B verify_followup.py <repo-root>
Exact integers / Fractions only. Reads sources; executes none of them.

Claims checked here (all conditional on the retained interfaces named in the citations):
 (A) The old ordinary leaf A0 (a0 = 384599/10^10, PR36 copied-centres bit supplier) has
     m=575, largest child 529, wire bits 28 (W < 2^28), least halving degree 9, row coefficient 252.
 (B) For its depth recurrence D(e)=0 (e<m), D(e)=1+D(529*floor(e/575)), we have
     D(e) <= 9*ceil(log2 max(2,e)) for ALL e >= 1 (PR141's 'fixed constant' c may be taken as 0;
     in fact D(e) <= 9*ceil(log2 e) - 81 for e >= 575).
     Hence log2 R0 < 28*D <= 252*ceil(log2 max(2,e)) and log_q R0 < (63/20)(1+t) for q > 2^80.
 (C) Stage envelopes sigma_j under three readings, compared with H_j for all 14 recorded stages.
"""
import hashlib, json, sys
from fractions import Fraction as Q
from pathlib import Path

R = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[2]


class Fail(Exception):
    pass


def need(c, m):
    if not c:
        raise Fail(m)


def line(path, n, text):
    ls = (R / path).read_bytes().decode().split('\n')   # grep-style numbering; read_text would turn bare CR into newlines
    need(n <= len(ls) and text in ls[n - 1], f'{path}:{n} must contain {text!r}')
    return f'{path}:{n}'


RET = 'research/five-stage-source527-banks/inputs/complex/retained/'
PR141 = 'research/coordinated-frames-and-entrance-banks/references/pr141/'
PR185 = 'research/community-round8-audit/pr185-proof-reviewed.md'
PR387 = 'research/bootstrap-depth-audit/sources/pr387/PROOF.md'
ROWS = 'research/five-stage-source527-banks/proof/three-stage-cover-rows.tex'
SPLIT = 'research/bootstrap-cutoff-time-space-split/'
CITES = [
    # A0's own supplier record (PR36 copied-centres assembly, retained in the source527 bundle).
    (RET + 'copied-centers-assembly.tex', 83, 'Take the strict finite savings $a=a_b=384599/10^{10}$'),
    (RET + 'copied-centers-assembly.tex', 38, '(t_b,w_b)=(9,28)'),
    (RET + 'copied-centers-assembly.tex', 40, 'where $m_j^{t_j}>2r_j^{t_j}$ for largest children $r_b=529'),
    (RET + 'copied-centers-assembly.tex', 43, '$D_j(e)=t_j\\lceil\\log_2e\\rceil$'),
    (RET + 'copied-centers-assembly.tex', 44, '$9\\cdot28+20\\cdot30=852$'),
    (RET + 'endpoint-gauge-complex.tex', 187, '$D(e)=0$ for $e<m$ and $D(e)=1+D(t_*\\lfloor e/m\\rfloor)$ otherwise.'),
    # PR141: the old-leaf row and its unquantified constant; its validator's definitions.
    (PR141 + 'ORIGINAL-PROOF.md', 9, 'old ordinary saving 384599/10^10'),
    (PR141 + 'ORIGINAL-PROOF.md', 25, 'W_complex^D_complex * W_coarse^D_coarse * W_old^D_old'),
    (PR141 + 'ORIGINAL-PROOF.md', 29, 'd*ceil(log2(max(2,e))) plus a fixed constant'),
    (PR141 + 'ORIGINAL-PROOF.md', 35, '| Old ordinary leaf | 575 | 529 | 28 | 9 | 252 |'),
    (PR141 + 'balanced_stopped.py', 52, 'while m**degree <= 2 * child**degree:'),
    (PR141 + 'balanced_stopped.py', 68, 'degree, bits = halving_degree(m, child), W.bit_length()'),
    # Reserve readings.
    (ROWS, 5, 'A sufficient common divisor is $W(w)^{D(n)}$ times the inherited old-leaf'),
    (ROWS, 7, 'logarithmic reserve'),
    (PR185, 9, 'Its ordinary wrapper borrows and restores its required high-digit stock internally.'),
    (PR185, 38, 'its new ordinary leaves need no *external* extra stock'),
    (PR185, 44, '**252 is not asserted to be the new leaf\'s actual internal reserve.**'),
    (PR387, 7, 'performs required high-digit borrowing internally'),
    (PR387, 7, 'It is retained conditionally, not proved anew here.'),
    (PR387, 9, 'ordinary calls decrease j'),
    (PR387, 36, 'The retained external 252 allowance is unused slack'),
    ('upstream/build/sections/04-swap.tex', 282, 'Rows for a general array can be supplied by a few of its own high digits.'),
    ('upstream/build/sections/04-swap.tex', 308, '\\rho=\\left\\lceil\\frac{k\\log W}{2\\log q}\\right\\rceil.'),
    ('upstream/build/sections/04-swap.tex', 354, 'by each required $W^j$'),
    # Selected producer (stage-1 package) and guard.
    ('research/bootstrap-depth-audit/sources/pr387/refine.py', 66, 'degraded_bit_saving=str(saving),bit_coarse=str(c)'),
    ('research/five-stage-source527-banks/proof/three-stage-cover-bit.tex', 198, 'Take $q>2^{80}$'),
    ('research/five-stage-source527-banks/proof/three-stage-cover-bit.tex', 212, 'The full role count is $W(w)=W_0q^{m^2(w-1)}$'),
    (SPLIT + 'README.md', 56, '(H3) `log_q S(e,w) ≤ σ·w·(1+log2 e)`'),
    (SPLIT + 'README.md', 94, 'valid exactly when every `σ_j ≤ 2^(H_j−3)`'),
    (SPLIT + 'verify.py', 48, 'A0 = Q(384599, 10**10)'),
    (SPLIT + 'verify.py', 182, "headroom = floor(L * row['gap'] / 2) - 4"),
]


def check_cites(cites):
    return [line(*c) for c in cites]


def old_leaf(m, r, bits, d, coef):
    """(A): least halving degree, coefficient, and W < 2^bits from bit_length semantics."""
    need(0 < r < m, 'child contracts')
    deg = 1
    while m ** deg <= 2 * r ** deg:
        deg += 1
    need(deg == d, f'least halving degree is {deg}, claimed {d}')
    need(coef == d * bits, f'row coefficient {coef} != {d}*{bits}')
    return deg


def depth(m, r, e):
    n = 0
    while e >= m:
        e, n = r * (e // m), n + 1
    return n


clog = lambda n: (max(2, n) - 1).bit_length()          # ceil(log2 max(2,n))


def depth_lemma(m, r, d, c, N):
    """(B): D(e) <= d*ceil(log2 max(2,e)) + c for e < N, plus the general proof's inequalities.

    Proof (all e): let e_0=e, e_{i+1}=r*floor(e_i/m), L=D(e) >= 1. Then e_i <= (r/m)^i e and
    m <= e_{L-1}. As (r/m)^d < 1/2, 2^floor((L-1)/d) <= e/m < 2^(ceil(log2 e) - 9) since m > 2^9.
    So floor((L-1)/d) <= ceil(log2 e) - 10 and L <= d*ceil(log2 e) - 9d. For e < m, D=0. Hence c=0.
    """
    need(m ** d > 2 * r ** d, 'd-level halving')
    need(m > 2 ** 9, 'm > 2^9 used in the tight form')
    need(all(depth(m, r, e) <= d * clog(e) + c for e in range(1, N)), f'depth bound with c={c}')
    need(all(depth(m, r, e) <= d * clog(e) - 9 * d for e in range(m, N)), 'tight form')
    return max(depth(m, r, e) - d * clog(e) for e in range(1, N))


def chain(bit, cx, a0):
    s = min(bit, (1 - Q(1, 10**24)) * cx) - Q(1, 10**27)
    a, rows = a0, []
    for _ in range(7):
        old, a = a, (1 - s) * s + s * a
        rows.append((a, min(s - a, 1 - a - s, 1 - a - s * (1 - old), 1 - s)))
    return rows


def fit_ok(L, g, ell):
    return L >= 2 and L * g * g >= 36 and L * g >= 2 * (4 + ell)


def headrooms(asm, a0):
    out = {}
    for key, side in (('bootstrap', 'complex_lower'), ('conservative_bootstrap', 'conservative_complex_lower')):
        rows, hs = chain(Q(asm['bit_lower']), Q(asm[side]), a0), []
        need(len(asm[key]) == 7, key + ' has 7 stages')
        for j, ((a, g), rec) in enumerate(zip(rows, asm[key]), 1):
            need(Q(rec['exponent']) == a and Q(rec['gap']) == g, f'{key}[{j}] chain from a0')
            L = int(rec['cutoff_log2'])
            H = (L * g / 2).__floor__() - 4
            need(fit_ok(L, g, H) and not fit_ok(L, g, H + 1), f'{key}[{j}] exact headroom')
            hs.append(H)
        out[key] = hs
    return out


def guard(sigma, H):
    """log2(8 sigma) <= H, sufficient check without forming 2^H: 8*num < 2^H <= 2^H*den."""
    x = 8 * Q(sigma)
    return x.numerator.bit_length() <= H


def main():
    cites = check_cites(CITES)
    # Selected-record bindings (as in the stage-1 package).
    asm_p = R / 'research/community-frontier1094-audit/ASSEMBLY.json'
    sel = json.loads((R / 'certificates/selected-result.json').read_text())
    need(hashlib.sha256(asm_p.read_bytes()).hexdigest() == sel['certificate_sha256'], 'selected binds ASSEMBLY')
    asm = json.loads(asm_p.read_text())
    exp_p = R / 'research/e8-bit-p111-tuned/tune/EXPECTED.json'
    need(asm['input_sha256']['package:inputs/PR589-EXPECTED.json'] == hashlib.sha256(exp_p.read_bytes()).hexdigest(), 'ASSEMBLY pins PR589')
    bit = json.loads(exp_p.read_text())['nine']['bit']
    need((bit['maxchild'], bit['halving_depth'], bit['row_coefficient']) == (42, 2, 18355), 'PR589 producer')
    need(18355 == 2 * (64 * 64 + 1) + 9909 + 252, 'PR589 external row = 8194 + 9909 + 252')
    # Structural twin of A0 with machine-readable fields (copied-reversed certificate, prior bit = a0).
    cr = json.loads((R / 'research/copied-reversed/certificate.json').read_text())
    need(cr['prior']['bit'] == '384599/10000000000', 'copied-reversed prior bit is a0')
    fb = cr['finite_bridge']['bit']
    need((fb['m'], fb['maxchild'], fb['wire_bits'], fb['halving_degree']) == (575, 529, 28, 9), 'm575 twin fields')
    need(2 ** 27 <= fb['W'] < 2 ** 28, 'twin W has 28 bits')

    # (A) + (B)
    A0 = Q(384599, 10**10)
    d = old_leaf(575, 529, 28, 9, 252)
    worst = depth_lemma(575, 529, d, 0, 1 << 21)
    # (B) consequence: log_q R0 < 28*D/80 <= 252(1+t)/80 = (63/20)(1+t) <= (63/20) w (1+t).
    RHO0 = Q(252, 80)
    need(RHO0 == Q(63, 20), 'rho0')

    # (C) per-stage envelopes. Parent per-stage term from the stage-1 package: 2*(64^2+1) = 8194.
    SIG = 2 * (64 * 64 + 1)
    need(SIG == 8194, 'parent coefficient')
    H = headrooms(asm, A0)
    need(H['bootstrap'][0] == H['conservative_bootstrap'][0] == 6562783, 'H_1')
    stages = []
    for key in ('bootstrap', 'conservative_bootstrap'):
        for j, h in enumerate(H[key], 1):
            s_int, s_alloc, s_hoist = Q(SIG), SIG + RHO0, SIG * j + RHO0
            for s in (s_int, s_alloc, s_hoist):
                need(guard(s, h), f'{key}[{j}] guard')
            stages.append(dict(record=key, stage=j, H=h, sigma_internal=str(s_int),
                               sigma_retained_252=str(s_alloc), sigma_hoisted=str(s_hoist),
                               log2_8sigma_hoisted_ceil=(8 * s_hoist).__ceil__().bit_length()))
    return dict(status='PASS_SELECTED_STOCK_ENVELOPE_ALL_RECORDED_STAGES_CONDITIONAL',
                old_leaf=dict(a0='384599/10^10', m=575, maxchild=529, wire_bits=28, halving_degree=d,
                              coefficient=252, c=0, worst_D_minus_bound_below_2_21=worst,
                              log_q_R0='< (63/20)(1+log2 e)'),
                sigma_internal=SIG, sigma_retained_252='163943/20', sigma_hoisted_j='8194*j + 63/20',
                stages=stages, citations=cites)


def controls():
    """Mutation controls: each must FAIL."""
    out = {}
    trials = {
        'wrong_halving_degree_8': lambda: old_leaf(575, 529, 28, 8, 224),
        'wrong_coefficient_251': lambda: old_leaf(575, 529, 28, 9, 251),
        'wrong_wire_bits_27': lambda: old_leaf(575, 529, 27, 9, 252),
        'negative_constant_c_minus_10': lambda: depth_lemma(575, 529, 9, -10, 1 << 16),
        'bigger_child_550_breaks_degree_9': lambda: depth_lemma(575, 550, 9, 0, 1 << 12),
        'wrong_citation_line': lambda: line(PR141 + 'ORIGINAL-PROOF.md', 34, '| Old ordinary leaf | 575 | 529 | 28 | 9 | 252 |'),
        'wrong_citation_text': lambda: line(RET + 'copied-centers-assembly.tex', 38, '(t_b,w_b)=(9,27)'),
        'wrong_a0_chain': lambda: headrooms(json.loads((R / 'research/community-frontier1094-audit/ASSEMBLY.json').read_text()), Q(384598, 10**10)),
        'sigma_above_headroom': lambda: need(guard(Q(2) ** (6562783 - 3) + 1, 6562783), 'guard'),
    }
    for name, f in trials.items():
        try:
            f()
            out[name] = 'UNEXPECTED_PASS'
        except Fail as e:
            out[name] = 'rejected: ' + str(e)[:90]
    need(all(v.startswith('rejected') for v in out.values()), 'all mutation controls rejected')
    return out


if __name__ == '__main__':
    try:
        res = main()
        res['controls'] = controls()
    except Fail as e:
        raise SystemExit('FAIL ' + str(e))
    print(json.dumps(res, indent=1))
