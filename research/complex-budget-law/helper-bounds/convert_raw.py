"""Convert a raw provider word (selected #1094 K8-SLOT-WORD; hpe/audit_local `run` semantics) into the
BORROW-WORD layout read by wordmodel.Word: no warps/removals, derived compensation (exact, as audit_local.run),
early-cleanup tail gates (cleanup[i] at frame U, run after B) as the 'restore' section.
usage: convert_raw.py RAW.json.gz OUT.json.gz"""
import sys, gzip, json
from fractions import Fraction as F
w = json.load(gzip.open(sys.argv[1]))
v, R = w['report']['v'], w['report']['R']; n = 2 * v + R


def apply(e, rows):
    if e[0] == 'add':
        t, s = e[1:3]; q = F(e[3], e[4]); dst = rows[t]
        for k, x in list(rows[s].items()):
            a = dst.get(k, F(0)) + q * x
            if a: dst[k] = a
            else: dst.pop(k, None)
    else:
        t = e[1]; q = F(e[2], e[3]); rows[t] = {k: q * x for k, x in rows[t].items()}


raw = [{k: F(1)} for k in range(n)]
for e in w['A'] + w['cut'] + w['B']:
    apply(e, raw)
dirt = [{} for _ in range(n)]
for y in range(v, 2 * v):
    dirt[y] = {k: a for k, a in raw[y].items() if k >= 2 * v}
for e in w['precondition']:
    apply(e, dirt)
comp = [['add', y, k, -a.numerator, a.denominator, []] for y in range(v, 2 * v) for k, a in sorted(dirt[y].items())]
tail = [w['cleanup'][i][:-1] + [U] for i, U in w.get('early_cleanup', [])]
out = dict(kind='RAW-AS-BORROW', ports=w['ports'], centers=w['centers'], center_frames=w['center_frames'],
           target_levels=w.get('target_levels', {}), start_extra=w.get('start_extra', {}),
           helper_final=w.get('helper_final', {}), removed=[], initial_corrections=[],
           precondition=w['precondition'], compensation=comp, A=w['A'], cut=w['cut'], B=w['B'],
           restore=tail, warps=[], report=w['report'], prelude_events=0)
json.dump(out, gzip.open(sys.argv[2], 'wt'))
print('compensation', len(comp), 'tail', len(tail), 'report', {k: w['report'][k] for k in ('local_profile', 'width_inventory') if k in w['report']})
