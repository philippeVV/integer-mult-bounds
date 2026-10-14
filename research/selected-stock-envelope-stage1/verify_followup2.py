#!/usr/bin/env python3
"""Selected stock envelope, follow-up 2: source-bound discharge of conditions (1), (3), (4).

Usage: python3 -I -B verify_followup2.py <repo-root>
Exact integers / Fractions only. Reads sources; executes none of them (every replay below is an
independent re-implementation, checked against the recorded values).

 (1) A0 = PR36's copied-centres bit supplier (a0 = 384599/10^10). Replay its physical two-stage
     histogram from the sha-pinned producer records + (25,23) corner certificate:
     re-derive the corner pivots/runs/data profile and all 315 zero-minor ranks, the copied-centre
     histogram replacement, and the full child histogram. Recompute m, max child, W, wire bits,
     least halving degree; compare with copied-centers-input/network, PR141's table and the
     retained assembly text (byte-identical to the sha-pinned notes/copied-centers-assembly.tex).
 (3) Role inventory of the selected PR589 bit word: exact replay of the tuned word (nested frames,
     F2 endpoint on every column) -> local paid profile; live = 2Kv T, physical = live + |ret|;
     nine-stage profile = 288*local + REM; rank-mass identity. The per-level role multiplier of
     W(w) is exactly physical = 28425 < 2^15, replacing the 52-bit overcount.
     Also re-checks the m-dependent finite inequalities of the GL_m(R_w) cover at m = 64.
 (4) Finite content of the fit lemma's (H1)/(H2) at all 14 stages: exact gaps/ordering, the H2
     exponent identity, mu_1 < 1 and mu_tau < 1 (tau = 1-s, with the 10^-16*32m^2 fallback envelope)
     for the selected m=64 producer, and mu < 1 at a0 for A0 (H1 base); H1 for j >= 2 follows from
     the lemma's conclusion (2) at j-1 given H2(j-1) and fit(j-1).
"""
import gzip, hashlib, json, sys
from collections import Counter
from fractions import Fraction as Q
from math import comb
from pathlib import Path

R = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[2]


class Fail(Exception):
    pass


def need(c, m):
    if not c:
        raise Fail(m)


def line(path, n, text):
    ls = (R / path).read_bytes().decode().split('\n')
    need(n <= len(ls) and text in ls[n - 1], f'{path}:{n} must contain {text!r}')
    return f'{path}:{n}'


def J(p):
    return json.loads((R / p).read_text())


def sha(p):
    return hashlib.sha256((R / p).read_bytes()).hexdigest()


RET = 'research/five-stage-source527-banks/inputs/complex/retained/'
PR141 = 'research/coordinated-frames-and-entrance-banks/references/pr141/'
SPLIT = 'research/bootstrap-cutoff-time-space-split/'
P589 = 'research/e8-bit-p111-tuned/'
GLOBALK = P589 + 'vendor/pr510/e8-r480-nine-stage/vendor/e8-r715-eleven-stage-whole/code/global_k.py'
MOMENT = P589 + 'vendor/pr510/e8-r480-nine-stage/vendor/e8-r740-replica32/code/engines/moment.py'
COVER = 'research/five-stage-source527-banks/proof/three-stage-cover-bit.tex'
FBA = 'research/five-stage-source527-banks/proof/FINITE_BIT_ACCOUNTING.md'
CCN = 'scripts/copied_centers_network.py'
CITES = [
    # (1) PR36 word -> parameters.
    (CCN, 26, 'AB = Q(384599, 10**10)'),
    (CCN, 41, 'W = 2*N+B1+B2'),
    (CCN, 59, 'z[m-2*h] += count'),
    (CCN, 83, 'z[1] += N  # Paid endpoint transform on the copied first output.'),
    (CCN, 89, 'maxchild=max(z)'),
    (CCN, 135, "db, dc = halving(bit['m'], bit['maxchild'])"),
    (CCN, 136, "wb, wc = bit['W'].bit_length()"),
    ('scripts/copied_centers/physical.py', 22, 'changed[1] += h'),
    ('scripts/copied_centers/physical.py', 23, 'changed[h] -= h'),
    ('scripts/copied_centers/corners.py', 105, 'blocks=[len(run) for run in runs if len(run)>1]+[m-2*d]'),
    ('notes/copied-centers-bit.tex', 14, '25&2300&45963&6925&1589&51299&600'),
    ('notes/copied-centers-bit.tex', 15, '23&1771&34764&5336&1324&38776&506'),
    ('notes/copied-centers-bit.tex', 18, 'R_h=c_h+q_h-|\\mathcal M_h|'),
    (RET + 'copied-centers-assembly.tex', 38, '(t_b,w_b)=(9,28)'),
    (RET + 'copied-centers-assembly.tex', 40, 'r_b=529'),
    (RET + 'copied-centers-assembly.tex', 83, 'a=a_b=384599/10^{10}'),
    (PR141 + 'ORIGINAL-PROOF.md', 35, '| Old ordinary leaf | 575 | 529 | 28 | 9 | 252 |'),
    # (3) PR589 role inventory and the cover.
    (GLOBALK, 47, 'K=(L-1)//2; T=a.replicas; M=1+7*L'),
    (GLOBALK, 51, 'DATA=2*K*V*T'),
    (GLOBALK, 72, 'LIVE=DATA'),
    (GLOBALK, 239, 'private_centre_reserve=9'),
    (GLOBALK, 248, 'physical=LIVE+9'),
    (P589 + 'tune/fastc.py', 12, 'REM={2:13920,4:14208,6:14208,7:21312,42:288}'),
    (P589 + 'tune/fastc.py', 13, 'W=28416;m=64;T=288'),
    (P589 + 'tune/EXPECTED.json', 119, '"physical": 28425'),
    (P589 + 'tune/EXPECTED.json', 123, '"bank_assignments": 0'),
    (COVER, 163, 'N_0=2m'),
    (COVER, 190, '6N_0(N_0-1)+3N_0+6(N_0-1)<32m^2.'),
    (COVER, 198, 'Take $q>2^{80}$'),
    (COVER, 211, 'W_0=|\\operatorname{GL}_m(\\mathbb F_q)|(2v+3R_{\\rm loc})'),
    (COVER, 212, 'W(w)=W_0q^{m^2(w-1)}'),
    (COVER, 225, '2m^3/q<10^{-16}'),
    # (4) fit lemma hypotheses and their finite ingredients.
    (SPLIT + 'README.md', 44, '(H1) `A_(j−1)` is a completed ordinary algorithm at every'),
    (SPLIT + 'README.md', 45, '(H2) Whenever `k(e) < e/2`, the stopped construction is a'),
    (SPLIT + 'README.md', 48, 'T(e)/V ≤ K·[e^(1−s) + e^(1−a\') + e^s(1+log2 e) + w^(1−a)]'),
    (SPLIT + 'README.md', 51, 'the weighted-induction factor `1+1/(1−μ_1)+1/(1−μ_τ)`'),
    (SPLIT + 'README.md', 64, '`T_j(e) ≤ C_j·V·e^(1−a\')` for all `e ≥ 1`'),
    (FBA, 79, 'adding the entire 10^-16*32m^2 fallback envelope per edge'),
    (FBA, 86, 'take A=C*(1+1/(1-mu_1)+1/(1-mu_tau))'),
    (FBA, 101, 'Starting with completed ordinary a0=384599/10^10'),
    ('research/bootstrap-depth-audit/sources/pr387/PROOF.md', 13, '`T_j(e)/V = O(e^(1-s) + e^((1-s)^2+s(1-a_(j-1))) + e^s log(e))`'),
    (MOMENT, 60, 'p=F(r*n,m*W)'),
    (MOMENT, 64, 'p=F(32*m*sum(H.values()),10**16*W)'),
]

# ---------------------------------------------------------------- exact real enclosures
GRID = 1 << 256
dn = lambda x: Q(x.numerator * GRID // x.denominator, GRID)
up = lambda x: -dn(-x)


def _atanh_ln(y, K=60):
    """ln y for 1 <= y < 2 via 2*atanh(z), z=(y-1)/(y+1) <= 1/3: (lower, upper)."""
    z = (y - 1) / (y + 1)
    s, t = Q(0), z
    for k in range(K):
        s += t / (2 * k + 1)
        t *= z * z
    return 2 * s, 2 * s + 2 * t / ((2 * K + 1) * (1 - z * z))


_l2 = _atanh_ln(Q(2))                     # z = 1/3
LN2 = (dn(_l2[0]), up(_l2[1]))


def ln_bounds(x):
    need(x >= 1, 'ln domain')
    k = 0
    while x >= 2:
        x /= 2
        k += 1
    lo, hi = _atanh_ln(x)
    return dn(lo + k * LN2[0]), up(hi + k * LN2[1])


def exp_bounds(lo, hi, n=30):
    need(0 <= lo <= hi < 1, 'exp domain')
    def poly(u):
        s = t = Q(1)
        for j in range(1, n + 1):
            t = t * u / j
            s += t
        return s, t
    a, _ = poly(lo)
    b, last = poly(hi)
    return dn(a), up(b + last * hi / (n + 1) / (1 - hi / (n + 2)))


def moment(H, m, W, a, bad=False):
    """(lower, upper) of sum_r (r n/(m W)) (m/r)^a [+ 32 m sum(H)/(10^16 W) m^a]  (moment.py:60,64 shape)."""
    lo = hi = Q(0)
    terms = [(r, n) for r, n in H.items()]
    if bad:
        terms.append(('bad', None))
    for r, n in terms:
        if r == 'bad':
            p, x = Q(32 * m * sum(H.values()), 10**16 * W), Q(m)
        else:
            need(0 < r < m and n > 0, 'child domain')
            p, x = Q(r * n, m * W), Q(m, r)
        l, u = ln_bounds(x)
        e, f = exp_bounds(a * l, a * u)
        lo += p * e
        hi += p * f
    return lo, hi


# ---------------------------------------------------------------- (1) PR36 word replay
def q_rank(rows):
    M = [list(map(Q, r)) for r in rows]
    rk, ncol = 0, len(M[0]) if M else 0
    for c in range(ncol):
        p = next((i for i in range(rk, len(M)) if M[i][c]), None)
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = M[i][c] / M[rk][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk


def corner_profile(cc, zero_minors=True):
    a, b = cc['dimensions']
    m, d = a * b, a + b - 1
    rows, cols = cc['R'], cc['C']
    need(len(rows) == len(cols) == d and cc['m'] == m and cc['d'] == d, 'corner shape')
    need(rows == (list(range(a)) + list(range(1, a)))[:d] and cols == ([a - 1] + list(range(a - 1)) + list(range(a)))[-d:], 'controlled boundary labels')
    left = [Q((r + 1) ** 2, sum(s * s for s in range(1, a + 1))) for r in range(a)]
    right = [Q((t + 1) ** 3, sum(s ** 3 for s in range(1, b + 1))) for t in range(b)]
    need([Q(x) for x in cc['witness']['left_weights']] == left and [Q(x) for x in cc['witness']['right_weights']] == right, 'witness weights')
    A = [[Q(int(rows[i] == cols[j])) / left[rows[i]] + Q(int(i % b == (m - d + j) % b)) / right[i % b] - 1
          for j in range(d)] for i in range(d)]
    free, piv, pvals = list(range(d)), [], []
    for i in range(d):
        j = max((c for c in free if A[i][c]), default=None)
        need(j is not None, 'singular corner')
        piv.append((i, j))
        pvals.append(A[i][j])
        free.remove(j)
        for k in range(i + 1, d):
            if A[k][j]:
                f = A[k][j] / A[i][j]
                A[k] = [x - f * y for x, y in zip(A[k], A[i])]
    need([(i, j) for i, j in piv] == [(w[0], w[1]) for w in cc['witness']['pivots']], 'recorded pivot positions')
    need(all(Q(w[2]) == pv for w, pv in zip(cc['witness']['pivots'], pvals)), 'recorded pivot values')
    runs = []
    for i, j in piv:
        if runs and runs[-1][-1] == [i - 1, j - 1]:
            runs[-1].append([i, j])
        else:
            runs.append([[i, j]])
    need(runs == cc['runs'], 'pivot runs')
    prof = dict(singletons=sum(len(r) == 1 for r in runs), blocks=[len(r) for r in runs if len(r) > 1] + [m - 2 * d], rank=m - d)
    if zero_minors:
        const, full = 1 << (a + b), (1 << (a + b + 1)) - 1
        def inc(edges, mask):
            bits = [1 << i for i in range(a + b + 1) if mask >> i & 1]
            return [[int(x == const or x == 1 << u or x == 1 << v) for x in bits] for u, v in edges] or [[0]]
        recs, pos, sr, sc = cc['zero_minor_certificates'], 0, [], []
        for i, pj in piv:
            for j in range(pj + 1, d):
                if j in sc:
                    continue
                rec = recs[pos]
                pos += 1
                need((rec['row'], rec['col'], rec['size']) == (i, j, i + 1), 'zero-minor coverage')
                re_ = [(rows[k], a + k % b) for k in sr + [i]]
                ce_ = [(cols[k], a + (m - d + k) % b) for k in sc + [j]]
                r1 = q_rank(inc(re_, rec['partition'])) if rec['partition'] else 0
                r2 = q_rank(inc(ce_, full ^ rec['partition'])) if full ^ rec['partition'] else 0
                need((r1, r2) == (rec['row_rank'], rec['column_rank']) and r1 + r2 < i + 1, 'zero minor forced')
            sr.append(i)
            sc.append(pj)
        need(pos == len(recs) == 315, 'all 315 zero minors')
    return prof


def producer_rows(axes, psi):
    out = []
    for row in axes:
        h = row['h']
        src = psi['producers'][str(h)]
        need((src['v'], src['additions'], src['outputs'], src['matches'], src['roles'], src['loss'], src['histogram']) ==
             (row['v'], row['c'], row['q'], row['matching'], row['R'], row['loss'], row['histogram']), f'h={h} record = partial-swap provenance')
        need(row['v'] == comb(h, 3) and row['R'] == row['c'] + row['q'] - row['matching'] and row['loss'] == h * (h - 1), f'h={h} counts')
        H = row['histogram']
        need(len(H) == h + 1 and sum(r * n for r, n in enumerate(H)) == h * row['R'] + 2 * row['loss'] and H[h] >= h, f'h={h} rank identity')
        Hc = list(H)
        Hc[1] += h
        Hc[h] -= h
        need(sum(r * n for r, n in enumerate(Hc)) == h * row['R'] + row['loss'], f'h={h} copied rank mass')
        out.append(dict(h=h, v=row['v'], R=row['R'], loss=row['loss'], H=Hc))
    return out


def bit_profile(rows, cp):
    (a, va, Ra, la, Ha), (b, vb, Rb, lb, Hb) = [(r['h'], r['v'], r['R'], r['loss'], r['H']) for r in rows]
    m, N = a * b, va * vb
    B = [N // va * Ra, N // vb * Rb]
    W, L = 2 * N + B[0] + B[1], N // va * la + N // vb * lb
    z = Counter()
    for h, cnt in ((a, B[0]), (b, B[1])):          # stage exteriors: one rank-h and one rank-(m-2h) child
        z[h] += cnt
        z[m - 2 * h] += cnt
    z[1] += cp['singletons'] * 2 * N
    for t in cp['blocks']:
        z[t] += 2 * N
    def ordinary(h, r, n):                        # ordinary interchange of a rank-r residual on h lines
        if 2 * r > h:
            z[2 * r - h] += n
            z[1] += (h - r) * n
        else:
            z[1] += r * n
    for h, v, H in ((a, va, Ha), (b, vb, Hb)):
        for r, n in enumerate(H):
            ordinary(h, r, N // v * n)
        ordinary(h, h - 1, 2 * N)
    z[1] += N
    z = {t: n for t, n in sorted(z.items()) if t and n}
    s = W * m - N + L
    need(all(0 < t < m for t in z) and sum(t * n for t, n in z.items()) == s, 'two-stage rank')
    return dict(m=m, N=N, B1=B[0], B2=B[1], W=W, L=L, total_rank=s, deficit=N - L, maxchild=max(z), child_multiplicities=z)


def halving_degree(m, r):
    d = 1
    while m ** d <= 2 * r ** d:
        d += 1
    return d


def part1(tamper=None):
    net = J('certificates/copied-centers-network.json')
    for p in ('certificates/copied-centers-bit-axes.json', 'certificates/copied-centers-corner-25-23.json',
              'certificates/copied-centers-input.json', 'notes/copied-centers-assembly.tex', 'notes/copied-centers-bit.tex',
              'scripts/copied_centers_network.py', 'scripts/copied_centers/physical.py', 'scripts/copied_centers/corners.py'):
        need(net['source_sha256'][p] == sha(p), 'PR36 network pins ' + p)
    need(sha(RET + 'copied-centers-assembly.tex') == net['source_sha256']['notes/copied-centers-assembly.tex'], 'retained tex = pinned PR36 note')
    need(net['assembly']['parameters']['a_bit'] == '384599/10000000000', 'PR36 a_bit = a0')
    axes, cc = J('certificates/copied-centers-bit-axes.json'), J('certificates/copied-centers-corner-25-23.json')
    if tamper:
        tamper(axes, cc)
    cp = corner_profile(cc)
    need(cp == cc['data_profile'], 'corner data profile')
    rows = producer_rows(axes, J('certificates/partial-swap-input.json'))
    p = bit_profile(rows, cp)
    inp = J('certificates/copied-centers-input.json')['bit']
    for k in ('m', 'N', 'B1', 'B2', 'W', 'L', 'total_rank', 'deficit', 'maxchild'):
        need(p[k] == inp[k] == net['bit']['counts'][k], 'PR36 replay ' + k)
    need({str(t): n for t, n in p['child_multiplicities'].items()} == inp['child_multiplicities'] == net['bit']['counts']['child_multiplicities'], 'child histogram')
    m, r, W = p['m'], p['maxchild'], p['W']
    d, bits = halving_degree(m, r), W.bit_length()
    fb = net['finite_bridge']['bit']
    need((m, r, W, d, bits) == (fb['m'], fb['maxchild'], fb['W'], fb['halving_degree'], fb['wire_bits']) == (575, 529, 188181929, 9, 28), 'A0 parameters')
    need(d * bits == 252, 'row coefficient 252')
    # H1 base (finite core): A0 contracts and its moment is strict at tau0 = 1 - a0.
    a0 = Q(384599, 10**10)
    mu1 = Q(p['total_rank'], W * m)
    lo, hi = moment(p['child_multiplicities'], m, W, a0)
    need(mu1 < 1 and hi < 1, 'A0 strict moments')
    return dict(m=m, maxchild=r, W=W, wire_bits=bits, halving_degree=d, row_coefficient=d * bits,
                N=p['N'], B1=p['B1'], B2=p['B2'], deficit=p['deficit'], corner_profile=cp,
                one_minus_mu1=str(1 - mu1), one_minus_mu_tau0_lower=str(1 - hi), one_minus_mu_tau0_lower_decimal=float(1 - hi))


# ---------------------------------------------------------------- (3) PR589 word replay and inventory
def frame_dims(w):
    return [q_rank([w['lines'][i] for i in f]) if f else 0 for f in w['frames']]


def replay589(w):
    V, N = w['v'], 2 * w['v'] + w['R']
    dim = frame_dims(w)
    mask = [sum(1 << i for i in f) for f in w['frames']]
    val, cur, H = [1 << i for i in range(N)], list(w['start']), Counter()
    need(len(cur) == len(w['final']) == N, 'register count 2v+R')
    def move(reg, f):
        a = cur[reg]
        if a != f:
            need(mask[a] & ~mask[f] == 0, f'frame nesting reg {reg}')
            if dim[f] > dim[a]:
                H[dim[f] - dim[a]] += 1
            cur[reg] = f
    def run(sec):
        for f, t, src in w[sec]:
            need(t not in src and len(set(src)) == len(src), 'distinct operands')
            for reg in (t, *src):
                move(reg, f)
            for s_ in src:
                val[t] ^= val[s_]
    run('D'); run('A')
    tops = [dim[cur[q]] for q in w['ret']]
    snap = [val[q] for q in w['ret']]
    for t, ks in enumerate(w['scatter']):
        for k in ks:
            val[V + t] ^= snap[k]
    run('B'); run('C')
    for reg, f in enumerate(w['final']):
        move(reg, f)
    for reg in range(N):
        need(val[reg] == ((1 << reg) ^ (1 << (reg - V)) if V <= reg < 2 * V else 1 << reg), f'F2 endpoint col {reg}')
    return dict(sorted((H + Counter(tops)).items()))


def part3(word_tamper=None):
    asm_p = 'research/community-frontier1094-audit/ASSEMBLY.json'
    need(J('certificates/selected-result.json')['certificate_sha256'] == sha(asm_p), 'selected binds ASSEMBLY')
    asm = J(asm_p)
    need(asm['input_sha256']['package:inputs/PR589-EXPECTED.json'] == sha(P589 + 'tune/EXPECTED.json'), 'ASSEMBLY pins PR589')
    E = J(P589 + 'tune/EXPECTED.json')
    wp = P589 + 'tune/bitE8mut-p111-tuned.json.gz'
    need(sha(wp) == E['word']['file_sha256'], 'tuned word file hash')
    w = json.loads(gzip.decompress((R / wp).read_bytes()))
    need(hashlib.sha256(json.dumps(w, sort_keys=True, separators=(',', ':')).encode()).hexdigest() == E['word']['content_sha256'], 'word content hash')
    if word_tamper:
        word_tamper(w)
    local = replay589(w)
    need({str(k): n for k, n in local.items()} == E['word']['local_histogram'], 'local paid profile')
    v, Rr, ret = w['v'], w['R'], w['ret']
    need((v, Rr, len(ret)) == (E['word']['v'], E['word']['R'], E['word']['retained_totals']) == (111, 439, 9), 'v, R, retained totals')
    need(sum(r * n for r, n in local.items()) == E['word']['legality']['rank_mass'], 'local rank mass')
    L = 9; K, T, m = (L - 1) // 2, 32, 1 + 7 * L
    need(m == 64, 'm = 1 + 7L')
    live = 2 * K * v * T                                  # global_k.py:51,72 with no bank stock
    physical = live + len(ret)                           # global_k.py:239,248 private centre reserve
    bit = E['nine']['bit']
    need(bit['bank_assignments'] == 0 and (live, physical) == (bit['live'], bit['physical']) == (28416, 28425), 'role inventory')
    REM = {2: 13920, 4: 14208, 6: 14208, 7: 21312, 42: 288}      # fastc.py:12 (retained hosted audit)
    P = Counter({r: T * L * n for r, n in local.items()})
    P.update(REM)
    P = dict(sorted(P.items()))
    need({str(k): n for k, n in P.items()} == E['nine']['profile'], 'nine-stage profile = 288*local + REM')
    mass = sum(r * n for r, n in P.items())
    need(mass == bit['rank_mass'] and max(P) == bit['maxchild'] == 42, 'rank mass / maxchild')
    need(m * live - mass == T * (2 * K * v - L * E['word']['decoder_constant']) == 10272, 'deficit identity')
    need(2 * 42 ** 2 < m ** 2 and halving_degree(m, 42) == bit['halving_depth'] == 2, 'two-level halving')
    # The cover at m = 64: log_q W(w) = log_q|GL_m(F_q)| + log_q(physical) + m^2(w-1) < m^2 w + log_q(physical).
    need(physical < 2 ** 15 and physical.bit_length() == 15, 'role multiplier has 15 bits')
    need(Q(15, 80) < 1, 'log_q multiplier < 1 at q > 2^80')
    overcount = physical * (1 + 32 * m * m * sum(P.values()))
    N0 = 2 * m
    need(6 * N0 * (N0 - 1) + 3 * N0 + 6 * (N0 - 1) < 32 * m * m, 'fallback count < 32 m^2 at m=64')
    need(Q(2 * m ** 3, 2 ** 80) < Q(1, 10**16), 'bad fraction 2m^3/q < 10^-16 at m=64')
    return dict(asm=asm, P=P, m=m, live=live, out=dict(
        word=dict(v=v, R=Rr, retained_totals=len(ret), local_paid_profile=local),
        live=live, physical=physical, role_multiplier=physical, role_multiplier_bits=15,
        previous_overcount_bits=overcount.bit_length(), nine_profile=P, rank_mass=mass, maxchild=42,
        log_q_W_w='< 4096 w + 15/80  (q > 2^80)',
        cover_inequalities_m64=dict(fallback=6 * N0 * (N0 - 1) + 3 * N0 + 6 * (N0 - 1), bound_32m2=32 * m * m,
                                    bad_fraction_num_over_2_80=2 * m ** 3)))


# ---------------------------------------------------------------- (4) fit-lemma hypotheses
def chain(bit, cx, a0):
    s = min(bit, (1 - Q(1, 10**24)) * cx) - Q(1, 10**27)
    a, rows = a0, []
    for _ in range(7):
        old, a = a, (1 - s) * s + s * a
        rows.append((old, a, min(s - a, 1 - a - s, 1 - a - s * (1 - old), 1 - s)))
    return s, rows


def part4(asm, P, m, live, a1, sigma_hoist=lambda j: 8194 * j + Q(63, 20)):
    a0 = Q(384599, 10**10)
    bl = Q(asm['bit_lower'])
    need(bl == Q(J(P589 + 'tune/EXPECTED.json')['nine']['bit_lower']), 'ASSEMBLY bit_lower = PR589 nine')
    mass = sum(r * n for r, n in P.items())
    mu1 = Q(mass, m * live) + Q(32 * m * sum(P.values()), 10**16 * live)      # a = 0: m^0 = 1
    need(mu1 < 1, 'selected mu_1 < 1 with fallback envelope')
    out = dict(one_minus_mu1_selected=str(1 - mu1), one_minus_mu1_decimal=float(1 - mu1), branches={})
    h1_prev = a1                                         # H1(1): A0 completed with strict moments (part 1)
    for key, side in (('bootstrap', 'complex_lower'), ('conservative_bootstrap', 'conservative_complex_lower')):
        s, rows = chain(bl, Q(asm[side]), a0)
        need(0 < a0 < s < bl, 'a0 < s < bit_lower')
        lo, hi = moment(P, m, live, s, bad=True)          # mu at tau = 1 - s
        need(hi < 1, key + ' mu_tau < 1')
        stages, h1 = [], h1_prev
        for j, ((aprev, a, g), rec) in enumerate(zip(rows, asm[key]), 1):
            need(Q(rec['exponent']) == a and Q(rec['gap']) == g, f'{key}[{j}] chain')
            need(aprev < a < s < 1 - a and g > 0, f'{key}[{j}] ordering/gaps')
            need((1 - s) ** 2 + s * (1 - aprev) == 1 - a, f'{key}[{j}] H2 exponent identity')
            need(s > a and 1 - a - s >= g and s * (1 - aprev) <= 1 - a - g, f'{key}[{j}] H2 term absorption')
            L = int(rec['cutoff_log2'])
            Hj = (L * g / 2).__floor__() - 4
            fit = L >= 2 and L * g * g >= 36 and L * g >= 2 * (4 + Hj) and (8 * Q(sigma_hoist(j))).numerator.bit_length() <= Hj
            need(fit, f'{key}[{j}] fit with hoisted sigma')
            h2 = True                                    # finite content: mu_1, mu_tau, gaps, identity (all checked)
            need(h1 and h2, f'{key}[{j}] H1/H2 finite content')
            stages.append(dict(stage=j, H1_finite=h1, H2_finite=h2, fit=fit, headroom=Hj))
            h1 = h1 and h2 and fit                       # lemma conclusion (2) at j gives H1 at j+1
        out['branches'][key] = dict(s=str(s), one_minus_mu_tau_lower=str(1 - hi), one_minus_mu_tau_decimal=float(1 - hi),
                                    weighted_induction_factor_upper=float(1 + 1 / (1 - mu1) + 1 / (1 - hi)), stages=stages)
    return out


def main():
    cites = [line(*c) for c in CITES]
    p1 = part1()
    p3 = part3()
    p4 = part4(p3['asm'], p3['P'], p3['m'], p3['live'], True)
    return dict(status='PASS_SELECTED_STOCK_FOLLOWUP2_CONDITIONS_1_3_4_FINITE',
                A0_replay=p1, PR589_inventory=p3['out'], fit_hypotheses=p4, citations=len(cites))


def controls():
    out = {}
    p3 = part3()
    asm, P = p3['asm'], p3['P']
    def bump_hist(ax, cc):
        ax[0]['histogram'][2] += 1; ax[0]['histogram'][1] -= 1
    def bad_label(ax, cc):
        cc['R'][5] = 6
    def zm_rank(ax, cc):
        cc['zero_minor_certificates'][100]['row_rank'] += 1
    def bad_axis(ax, cc):
        ax[1]['h'] = 22
    def frame_swap(w):
        w['A'][100][0] = w['A'][101][0] if w['A'][101][0] != w['A'][100][0] else 0
    def gate_doubled(w):
        w['A'].insert(500, list(w['A'][500]))
    trials = {
        'pr36_histogram_shift': lambda: part1(bump_hist),
        'pr36_corner_label': lambda: part1(bad_label),
        'pr36_wrong_axis_dim': lambda: part1(bad_axis),
        'pr36_zero_minor_rank': lambda: part1(zm_rank),
        'pr36_claim_W_2pow28_wire_bits_28': lambda: need((2 ** 28).bit_length() == 28, 'wire bits of 2^28 is 29'),
        'pr36_claim_halving_degree_8': lambda: need(halving_degree(575, 529) == 8, 'halving degree is 9'),
        'pr589_gate_frame_moved': lambda: part3(frame_swap),
        'pr589_gate_doubled': lambda: part3(gate_doubled),
        'physical_claim_28426': lambda: need(2 * 4 * 111 * 32 + 9 == 28426, 'physical is 28425'),
        'mu_tau_at_2x_bit_lower': lambda: need(moment(P, 64, 28416, 2 * Q(asm['bit_lower']), True)[0] < 1, 'mu lower bound >= 1'),
        'bad_fraction_q_2pow60': lambda: need(Q(2 * 64 ** 3, 2 ** 60) < Q(1, 10**16), 'bad fraction'),
        'sigma_hoisted_above_headroom': lambda: part4(asm, P, 64, 28416, True, lambda j: Q(2) ** 6562781),
        'H1_base_false': lambda: part4(asm, P, 64, 28416, False),
        'live_W_28000_breaks_mu1': lambda: part4(asm, P, 64, 28000, True),
        'wrong_citation_line': lambda: line(COVER, 211, '6N_0(N_0-1)+3N_0+6(N_0-1)<32m^2.'),
        'wrong_citation_text': lambda: line(GLOBALK, 248, 'physical=LIVE+8'),
    }
    for name, f in trials.items():
        try:
            f()
            out[name] = 'UNEXPECTED_PASS'
        except (Fail, IndexError, KeyError) as e:
            out[name] = 'rejected: ' + str(e)[:90]
    need(all(v.startswith('rejected') for v in out.values()), 'all mutation controls rejected: ' + json.dumps(out))
    return out


if __name__ == '__main__':
    try:
        res = main()
        res['controls'] = controls()
    except Fail as e:
        raise SystemExit('FAIL ' + str(e))
    print(json.dumps(res, indent=1, default=str))
