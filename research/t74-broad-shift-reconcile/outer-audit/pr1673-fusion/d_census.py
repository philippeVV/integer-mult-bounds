"""Part D: census of every cost/side term written in PR1673 (note + layout/tensor/support/halo notes), evaluated
exactly at PR2044's point (eps=1-2u, beta=1-6u, xi from its ledger, theta~p^-2, halo p^48, core p^128, W=p^(1/2),
Strassen tau=45/16). Each row: exponent e of p in cost/(V p^(1-kappa)) form, i.e. we print e_cost - (1-kappa);
a cost row is paid iff that is < 0. Side conditions print the quantity that must be < 0 (or > 0 as labelled).
Also the PR1673-regime values (halo 3/core 8, theta~1/d, PR1673's own eps/beta/xi) for comparison."""
import json, sys
from fractions import Fraction as F
from decimal import Decimal, getcontext
getcontext().prec = 60
sys.set_int_max_str_digits(0)
spec = json.load(open(sys.argv[1]))
DIG = spec['outer_refinement']['endpoint_digits']
u = F(1, 10**(DIG + 10))
b = F(spec['expected']['complex_lower']) - u; q = b - u; a = (q + b) / 2
def point(eps, beta, xi, halo, core, thp, kappa):
    return dict(eps=eps, beta=beta, xi=xi, halo=halo, core=core, thp=thp, kappa=kappa)
c = F(3, 32); tau = F(45, 16)
eps = 1 - 2 * u; beta = 1 - 6 * u; xi = eps * (1 + q) - c + u
kap = eps * q - F(1, 10**(DIG + 22))
P2044 = point(eps, beta, xi, 48, 128, 2, kap)
# PR1673 T74 Strassen branch (note sec. t74): q_T, r0=1e-6, zeta=1e-8, s=1e-5
qT = F('0.00167971337478') - F(1, 10**10); r0 = F(1, 10**6); ze = F(1, 10**8); s0 = F(1, 10**5)
e1 = (1 + qT + qT * (c - s0 - r0 - ze)) / (1 + qT)**2; x1 = e1 * (1 + qT) - c + s0; b1 = 1 - x1 - r0 - ze
P1673 = point(e1, b1, x1, 3, 8, None, e1 * qT - F(1, 10**11))
def rows(P, qq):
    e, bt, x, H, C, th, k = (P[z] for z in ('eps', 'beta', 'xi', 'halo', 'core', 'thp', 'kappa'))
    h = 1 - e                                   # log s_i = p^h
    cell = th if th is not None else e          # phase-cell length exponent (1/theta)
    w = F(1, 2)
    R = {}
    one = 1 - k
    def cost(name, ex, where): R[name] = ('cost', ex - one, where)
    def side(name, ex, where): R[name] = ('side<0', ex, where)
    cost('main transform V p d^-q', 1 - e * qq, 'note:290,298')
    cost('global routing (incl. one global low-bit gather)', 1 - a, 'note:299; OUTER-AUDIT:52')
    cost('fused Gaussian (d/g) g^(1-q) polylog H calls', e - qq * bt, 'note:136,300; TENSOR-COMPILER:147-151')
    cost('grouped bit routing (d/g) g^(1-a) per H call', e - a * bt, 'note:301; SUPPLIER-SUPPORT:86-87')
    cost('group scans: polylog diag scans, 2^g shifts, q-embeds, union gather/scatter, halo scan', e - bt, 'note:254,264,302')
    cost('boundary algebra (d/B) W^(omega-1)', e - x + (tau - 1) / 2, 'note:241,303')
    cost('PR1673 union, U scanned once per axis: (d/g) g^2 W/B', e + bt + w - x, 'note:250-252,304 (replaced in PR2044)')
    cost('union transposes/merge (d/g) g mu polylog', e + w - x, 'PHASE-LAYOUT-AUDIT:103-108; UNION-SCHEDULE:96')
    cost('core-face part of E_i: (d/g) g^2 W/L', e + bt + w - C, 'PHASE-LAYOUT-AUDIT:44-45')
    cost('dyadic-tail boundary W log B per cell, times W^(omega-2)', e + w - cell + (tau - 2) / 2, 'note:209-212')
    cost('packed polynomial products', 1 - e, 'note:306')
    cost('core-face corrections (d/g)(gA/L) g p^2', e + bt + H - C + 2, 'note:268-272,307')
    cost('halo sparse radix key passes (d/g)(gA/L) K0, K0=g h', e + H - C + bt + h, 'HALO-THEOREM:150-156')
    cost('halo metadata K0^2/B fallback (d/g)(gA/L)K0^2/p', e + H - C + 2 * (bt + h) - 1, 'HALO-THEOREM:41-42')
    side('membership g h^2 per p-bit record', bt + 2 * h - 1, 'note:254; PHASE-LAYOUT-AUDIT:123-130')
    side('union header g log L per record', bt - 1, 'PHASE-LAYOUT-AUDIT:119-121')
    side('group key g log s per record', bt + h - 1, 'note:273')
    side('tensor prefix guard g polylog bits vs p', bt - 1, 'TENSOR-COMPILER:21,110')
    side('H0 -> H rescale guard O(g) bits', bt - 1, 'TENSOR-COMPILER:154-157')
    side('support span gK, reservoirs O(g log^2 p) vs n=Theta(p)', bt - 1, 'SUPPLIER-SUPPORT:80,121')
    side('K-neighbourhood (log p) inside one axis (log s=p^h): -h', -h, 'note:126')
    side('core side log L inside one axis: -h', -h, 'OUTER-AUDIT:52')
    side('cap B below phase cell: xi - cell', x - cell, 'note:213 (B=o(d)); ledger cap_below_*')
    side('halo volume gA/L -> 0', bt + H - C, 'HALO-THEOREM:28,93-101')
    side('PR1673 precision guard g a B (linear gauge) vs p', bt + x - 1, 'note:276-281 (replaced in PR2044)')
    if th is not None:
        side('PR1758/1786 chirp guard g a theta B^2 vs p', bt + 2 * x - th - 1, 'PR2044 ledger group_gauge_power_slack')
    side('global scalar range d u vs p', e - 1, 'note:283')
    side('tables gN per coarse block vs N^g: (gN)/(N^g) [exponent of log]', F(-1), 'TENSOR-COMPILER:130 (superpolynomially small)')
    return R
spec2 = {}
for name, P, qq in [('PR2044', P2044, q), ('PR1673', P1673, qT)]:
    R = rows(P, qq)
    spec2[name] = {k: dict(kind=v[0], value=format(Decimal(v[1].numerator) / Decimal(v[1].denominator), '.6E'), ok=bool(v[1] < 0), where=v[2]) for k, v in R.items()}
for name in spec2:
    print('==', name)
    for k, v in spec2[name].items(): print(f"  {'OK ' if v['ok'] else 'BAD'} {v['kind']:6s} {v['value']:>16s}  {k}   [{v['where']}]")
json.dump(spec2, open(sys.argv[2], 'w'), indent=1)
