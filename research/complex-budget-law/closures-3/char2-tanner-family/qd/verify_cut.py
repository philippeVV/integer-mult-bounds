"""expQD: independent check of the balanced bisection of PG(2,16) - M (Singer model, point i ~ line i+d, d in D) with cut 1542,
and of the cut+viol consequence that keeps it from being the middle of an exact 772 schedule.
Lemma (expQD): for any schedule and any gap with prefix S (cut c), |R| >= (c + viol)/2, viol = sum_v [(e_part(v, earlier) - 8)^+ +
(e_part(v, later) - 8)^+] over the order induced on each side.  The first vertex of a side has all its internal neighbours later,
the last all of them earlier; if every vertex has internal degree >= 9, each side contributes >= 2, so viol >= 4.
An exact schedule (|R| = S = 772) through this gap would need c + viol <= 1544, but 1542 + 4 = 1546.
usage: verify_cut.py [file.json]   (default ../data/singsa_59_1.json, converted from expQD singsa_59_1.pkl)"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from singer import singer
n, D, _ = singer()
f = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '..', 'data', 'singsa_59_1.json')
d = json.load(open(f)); S = [int(x) for x in d['S']]; mate = [int(x) for x in d['mate']]
assert len(S) == 2 * n and sum(S) == n
M = {(i, mate[i] - n) for i in range(n)}
assert len(M) == n and len({j for _, j in M}) == n and all((j - i) % n in D for i, j in M)
E = [(i, (i + x) % n) for i in range(n) for x in D if (i, (i + x) % n) not in M]
assert len(E) == 16 * n
deg = [0] * (2 * n); intdeg = [0] * (2 * n)
for i, j in E:
    deg[i] += 1; deg[n + j] += 1
    if S[i] == S[n + j]: intdeg[i] += 1; intdeg[n + j] += 1
assert set(deg) == {16}
cut = sum(1 for i, j in E if S[i] != S[n + j])
mins = [min(intdeg[v] for v in range(2 * n) if S[v] == side) for side in (0, 1)]
print(f"{os.path.basename(f)}: balanced {sum(S)}/{2 * n - sum(S)}, 16-regular, perfect Singer-type matching, cut {cut}")
print(f"  min internal degree per side {mins} -> viol >= 4 -> cut + viol >= {cut + 4} > 1544 = 2S: not the middle of an exact 772 schedule")
ok = cut == 1542 and min(mins) >= 9
print('PASS_EXPQD_CUT' if ok else 'FAIL_EXPQD_CUT')
