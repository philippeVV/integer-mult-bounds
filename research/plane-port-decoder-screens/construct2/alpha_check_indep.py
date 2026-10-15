import sys, json, gzip, hashlib, tarfile
from ortools.sat.python import cp_model
TAR = sys.argv[1]; R = json.load(open(sys.argv[2]))['r']
with tarfile.open(TAR) as t:
    raw = t.extractfile('candidate/strict/global/base/providers/gauge44/CANONICAL-WORD.json.gz').read()
assert hashlib.sha256(raw).hexdigest().startswith('16bc0ebb04b4be6c')
P = json.loads(gzip.decompress(raw))['ports']
P = [p if isinstance(p, int) else int(''.join(map(str, p)), 2) for p in P]
d = lambda a, b: bin(a & b).count('1') & 1
assert len(P) == len(R) == 112 and len(set(P)) == 112
for q, r in zip(P, R): assert d(q, q) == 1 and d(r, q) == 0 and d(r, r) == 1
orth = lambda a, b: all(d(x, y) == 0 for x in (P[a], R[a]) for y in (P[b], R[b]))
E = [(a, b) for a in range(112) for b in range(a + 1, 112) if orth(a, b)]
m = cp_model.CpModel(); x = [m.NewBoolVar('') for _ in range(112)]
for a, b in E: m.AddBoolOr([x[a].Not(), x[b].Not()])
m.Maximize(sum(x)); s = cp_model.CpSolver(); s.parameters.num_workers = 2; s.parameters.max_time_in_seconds = 1800
st = s.Solve(m); print('edges', len(E), s.StatusName(st), 'alpha', s.ObjectiveValue(), 'bound', s.BestObjectiveBound())
