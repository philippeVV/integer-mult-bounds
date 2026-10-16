"""Per-q table from c1_all.log: best char-0 S over families with MDS local codes, best char-2 value of the
same frequency sets, the expQK conversion a = 1 - ln((E-S)/S)/ln 5 (girth-6, b = 5)."""
import json, math, collections
E = {}; best0 = {}; best2 = {}; q = None
N = lambda q: q*q + q + 1
for line in open(__file__.rsplit('/', 1)[0] + '/c1_all.log'):
    if line.startswith('== q='): q = int(line[5:]); continue
    if not line.startswith('{'): continue
    r = json.loads(line); n = q if q % 2 == 0 else q + 1; E[q] = N(q)*n
    mds = r['nonMDS_A'].startswith('0/') and r['nonMDS_B'].startswith('0/')
    if not mds: continue
    d = best0 if r['char'] != 2 else best2
    if r['S'] > d.get(q, (0, ''))[0]: d[q] = (r['S'], r['fam'])
a = lambda S, E: (1 - math.log((E - S)/S)/math.log(5)) if S else float('-inf')
print(f"{'q':>3} {'N':>4} {'E':>5} | {'best char-0 S (MDS)':>26} {'S/E':>6} {'a':>7} | {'best char-2 S (MDS)':>32} {'S/E':>6}")
for q in sorted(E):
    s0, f0 = best0.get(q, (0, '-')); s2, f2 = best2.get(q, (0, '-'))
    print(f"{q:>3} {N(q):>4} {E[q]:>5} | {s0:>5} {f0:>20} {s0/E[q]:6.4f} {a(s0,E[q]):7.3f} | {s2:>5} {f2:>26} {s2/E[q]:6.4f}")
ok = all(best0.get(q, (0,))[0] / E[q] < 1/6 and best2.get(q, (0,))[0] >= best0.get(q, (0,))[0] for q in E)
print('PASS_EXPRB_TABLE' if ok else 'FAIL_EXPRB_TABLE')   # every char-0 S/E < 1/6 (a < 0), and char 2 is never below char 0
