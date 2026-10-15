"""Fault controls: each mutated word must be rejected by the checker.  Usage: controls.py REPO WORD.json OUTDIR"""
import json, os, subprocess, sys, copy
from fractions import Fraction as Q
REPO, WORD, OUT = sys.argv[1:4]; w = json.load(open(WORD)); ops = w['ops']
chk = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'check_p104_relay_construction.py')
def idx(pred, k=0): return [i for i, o in enumerate(ops) if pred(o)][k]
muts = {}
m = copy.deepcopy(w); i = idx(lambda o: o[0] == 'shear', 5); m['ops'][i][3] = str(Q(m['ops'][i][3]) + 1); muts['shear_coefficient'] = m
m = copy.deepcopy(w); i = idx(lambda o: o[0] == 'read'); m['ops'][i][2] = (m['ops'][i][2] + 1) % 104; muts['read_wrong_source_frame'] = m
m = copy.deepcopy(w); i = idx(lambda o: o[0] == 'write' and o[4] == 0); del m['ops'][i]; muts['drop_dirt_compensation'] = m
m = copy.deepcopy(w); i = idx(lambda o: o[0] == 'move' and len(o[2]) == 7); m['ops'][i][2] = m['ops'][i][2][:6] + [511 ^ m['ops'][i][2][6]]; muts['non_nested_move'] = m
m = copy.deepcopy(w); i = idx(lambda o: o[0] == 'shear', -1); del m['ops'][i]; muts['drop_restoration_shear'] = m
m = copy.deepcopy(w); i = idx(lambda o: o[0] == 'write' and o[4] > 0); m['ops'][i][4] = 0; muts['early_target_write'] = m
ok = True
for name, mw in muts.items():
    p = os.path.join(OUT, 'ctl_%s.json' % name); json.dump(mw, open(p, 'w'))
    r = subprocess.run([sys.executable, '-I', '-B', chk, REPO, p], capture_output=True, text=True)
    rej = r.returncode != 0 and 'FAIL' in (r.stdout + r.stderr)
    print('%-26s %s  %s' % (name, 'REJECTED' if rej else 'ACCEPTED!!', (r.stdout + r.stderr).strip().splitlines()[-1][:90]))
    ok &= rej; os.remove(p)
print('CONTROLS_PASS' if ok else 'CONTROLS_FAIL')
