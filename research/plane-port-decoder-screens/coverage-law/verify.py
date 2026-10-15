"""expSB verifier: exact claims of README.md (coverage-law scouting). Prints PASS_EXPSB at the end.
Floating-point LM screens (c_lm.py, e_planes.py) are NOT re-run here; they are evidence only."""
import json
from fractions import Fraction as Q
from common import *
import f8
from ortools.sat.python import cp_model
P = load_ports(); n = len(P); base = cap(WORD)
def be(u, target=None):
    target = base if target is None else target
    lo, hi = 0, 600
    while hi - lo > 1:
        mid = (lo + hi) // 2
        try: ok = cap(WORD, centre=mid, unpaid=u) > target
        except Exception: ok = False
        lo, hi = (mid, hi) if ok else (lo, mid)
    return lo
# 1. coverage-law pricing: break-even and x1.3 thresholds, gain needed above the 2-integral floor tau
for tau, want in ((224, 163), (336, 263)):
    u = Q(tau, 112); b = be(u); assert b == want, (tau, b)
    print('tau', tau, 'break-even C', b, 'gain needed', tau - b, ' x1.3 needs C <=', be(u, Q(13, 10) * base), ' x2 needs C <=', be(u, 2 * base))
# 2. reductio of multi-total (matrix-coefficient) pricing: one query over all 112 line ports costs 9
assert pair_rank(P, P) == 9
loc = {r: c for r, c in WORD.items() if r != 8}; loc[9] = 1
r = cap(loc, centre=9) / base; assert r > Q(25, 10); print('line word, G/2 as one multi-total query of cost 9: ceiling x%.3f' % float(r))
# 3. paired identity fit: PR2126 invoice 112 vs rank-1-term invoice 224 (fails: no positive saving)
cert = json.load(open(MAIN + '/research/paired-plane-port-matching-seed6c563a/CERTIFICATE.json'))
inv = [pair_rank([P[i], P[j]], [P[i], P[j]]) for i, j in cert['matching']]
assert sum(inv) == 112 and all(v == 2 for v in inv)
try:
    cap(WORD, 224, 2); raise SystemExit('224 should fail')
except ValueError:
    pass
print('paired identity: matrix-term invoice 112, rank-1-term invoice 224 -> rejected by ceiling.py')
# 4. alpha of the plane-orthogonality graph: paired planes (48) and one F8 structure (40)
def alpha(free):
    m = cp_model.CpModel(); x = [m.NewBoolVar('') for _ in range(n)]
    for a in range(n):
        for b in range(a + 1, n):
            if free[a][b]: m.Add(x[a] + x[b] <= 1)
    m.Maximize(sum(x)); s = cp_model.CpSolver(); s.parameters.num_workers = 2; s.parameters.max_time_in_seconds = 300
    assert s.Solve(m) == cp_model.OPTIMAL; return int(s.ObjectiveValue())
pl = {}
for i, j in cert['matching']: pl[i] = pl[j] = (P[i], P[j])
fp = [[a != b and all(dot(x, y) == 0 for x in pl[a] for y in pl[b]) for b in range(n)] for a in range(n)]
ap = alpha(fp); assert ap == 48; print('paired planes: alpha(Gamma_E) =', ap, '(generic plane choices reach 23)')
F, ff = f8.build(1, 6)
assert all(f8.tr(f8.h(F[a], F[b])) == dot(P[a], P[b]) for a in range(n) for b in range(n)), 'trace form != dot'
af = alpha(ff); assert af == 40; print('F8-line ports (trial 6): free pairs', sum(map(sum, ff)) // 2, 'alpha =', af)
# 5. decoupled designs: I_S (+) restricted incumbent; cost per upgraded port >= dim E_a > (8/9)(dim E_a - 1)
for u in range(2, 10): assert u > Q(8, 9) * (u - 1)
print('PASS_EXPSB')
