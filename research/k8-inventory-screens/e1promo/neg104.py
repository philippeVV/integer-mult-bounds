"""Audit-level negative check of the ceiling claim: for each E1 helper whose exact ceiling is its entrance, build E2 words
(minimal completion, legality ignored) and require audit_local.run to reject them.  usage: neg104.py WORD cones OUT.json [per-helper]
(expB)"""
import sys, os, json, gzip, io, contextlib, importlib.util, random, tempfile
from multiprocessing import Pool
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
TMPROOT = os.path.join(os.path.dirname(os.path.abspath(sys.argv[3])), 'neg104-tmp'); os.makedirs(TMPROOT, exist_ok=True)
import e1lib as L
AUD = str(__import__('pathlib').Path(__file__).resolve().parents[3]) + '/research/gauge44-leaver8/independent/audit_local.py'
w, sha = L.load(sys.argv[1]); W = L.Word(w); W.sha = sha; W.prepare(); W.compute_cones(sys.argv[2]); C, gc = W.ceilings()
per = int(sys.argv[4]) if len(sys.argv) > 4 else 2
E1 = [q for q in W.entering if L.dim(W.start[q]) == 1 and C[q] == W.start[q]]
rng = random.Random(1148)
jobs = []
for q in E1:
    E = W.start[q]; opts = []
    seen = set()
    for x in range(1, 512):
        Ep = L.add(E, L.mask_of_basis((x,)))
        if Ep == E or Ep in seen: continue
        seen.add(Ep); raised, dprof, legal, bad, touched = W.complete({q: Ep}); assert not legal
        opts.append((len({b[0] for b in bad if b[0] != 'gate_ceiling'}), x, raised))
    opts.sort(key=lambda o: (o[0], o[1]))
    pick = [opts[0]] + rng.sample(opts[1:], per - 1)
    for _, x, raised in pick:
        jobs.append((q, x, {gi: list(L.basis_of_mask(m)) for gi, m in raised.items()}))

def run(job):
    q, x, raised = job
    s = importlib.util.spec_from_file_location('audit_local', AUD); A = importlib.util.module_from_spec(s); s.loader.exec_module(A)
    nw = json.loads(json.dumps(w))
    nw['start_extra'][str(q)] = list(L.basis(list(nw['start_extra'][str(q)]) + [x]))
    for gi, fr in raised.items():
        sec, i = W.ev[gi][0], W.ev[gi][1]
        if sec in ('A', 'B'): nw[sec][i][-1] = fr
        else: nw['early_cleanup'][i][1] = fr
    d = tempfile.mkdtemp(dir=TMPROOT); p = os.path.join(d, 'w.json.gz')
    open(p, 'wb').write(gzip.compress(json.dumps(nw, separators=(',', ':')).encode(), mtime=0))
    from pathlib import Path
    try:
        with contextlib.redirect_stdout(io.StringIO()): A.run(Path(p), Path(d) / 'r.json')
        got = 'ACCEPTED'
    except (RuntimeError, AssertionError) as e:
        got = str(e.args[0])[:60] if e.args else type(e).__name__
    return q, x, got
import multiprocessing as _mp
with _mp.get_context("fork").Pool(6) as pool: res = pool.map(run, jobs)
kinds = Counter(g.split(' (')[0] for _, _, g in res)
print('helpers', len(E1), 'words audited', len(res), 'rejections by reason', dict(kinds))
json.dump(dict(word=sha, helpers=E1, results=res), open(sys.argv[3], 'w'))
assert all(g != 'ACCEPTED' for _, _, g in res)
print('PASS every forced E2 word for a ceiling-1 helper is rejected by audit_local')
