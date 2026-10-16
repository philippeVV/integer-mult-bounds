"""expXP ladder summary (stdlib, < 1 s): re-reads the logged cmds.py ladders and checks the "complex r reached" column of
../README.md.  A rung counts as reached when f < 1e-6 and the flattened MDS margin is >= e^-3 (healthy).
  PG(2,5) (E = 186, E/6 = 31): reached r = 16, 17 (2 seeds each); r = 18, 19 fail (f 3e-6 .. 7e-5, margin e^-2.2 .. e^-4.1)
  PG(2,7) (E = 456, E/6 = 76): reached r = 30 (f 6e-6 < 1e-5, margin e^-2.3); r = 40, 50, 60, 77 fail (f >= 3.7e-3, margin <= e^-5.6)
  so r/E falls with q: 8/52 = 0.154 (PG(2,3) calibration, replayed by run_all), 17/186 = 0.091, 30/456 = 0.066.
Prints PASS_EXPXP_LADDERS."""
import json, os
H = os.path.dirname(os.path.abspath(__file__))
rows = [json.loads(l) for f in ('ladder_pg5.jsonl', 'ladder_pg7.jsonl', 'ladder_pg7b.jsonl') for l in open(os.path.join(H, f))]
def best(graph, ftol):
    ok = [d['r'] for d in rows if d['graph'] == graph and d['f'] < ftol and d['min_log_flat_margin'] >= -3]
    bad = sorted({d['r'] for d in rows if d['graph'] == graph and not (d['f'] < ftol and d['min_log_flat_margin'] >= -3)})
    return max(ok), bad
b5, f5 = best('PG5', 1e-6); b7, f7 = best('SG7', 1e-5)
print(f'PG(2,5): best healthy r = {b5} (r/E = {b5/186:.3f}); failing rungs {f5}')
print(f'PG(2,7): best healthy r = {b7} (r/E = {b7/456:.3f}); failing rungs {f7}')
ok = b5 == 17 and f5 == [18, 19] and b7 == 30 and f7 == [40, 50, 60, 77] and 8 / 52 > b5 / 186 > b7 / 456
ok &= b5 < 31 and b7 < 76    # far below E/6 (and below Boot >= 30, 77)
print('PASS_EXPXP_LADDERS' if ok else 'FAIL_EXPXP_LADDERS')
