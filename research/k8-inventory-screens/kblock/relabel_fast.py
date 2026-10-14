"""Fast relabelling sweep (cone model, see relabel_cone.py): all 8 k8 unit circuits are exact relabellings of cell 01+'s
unit (checked below), and the cone model is identical for all cells (premises checked by relabel_cone.check_premises),
so the 720 paired relabellings of the 01+ unit cover every cell's relabelling space.
For every sigma and every leaver source x_j: x_j can stay inside its block frame only if its unique x<-x gate is with its
block partner (otherwise the partner's cone carries another block's source); for the survivors, an exhaustive search
over the trace-equivalent orders (exact same linear map) decides lowness together with frame legality.
usage: python3 relabel_fast.py CODED-CELLS-k8.json"""
import sys, json, itertools, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from relabel_cone import search, idx

cells = json.loads(Path(sys.argv[1]).read_text())
base = [(t, s) for t, s, q in cells[0]['unit']]
rel = lambda u, sig: [(t[0] + str(sig[int(t[1])]), s[0] + str(sig[int(s[1])])) for t, s in u]
for c in cells:
    u = [(t, s) for t, s, q in c['unit']]
    assert any(rel(base, s) == u for s in itertools.permutations(range(6))), c['cell']
print('all 8 unit circuits are exact paired relabellings of the 01+ unit', flush=True)
t0 = time.time(); stats = dict(sigmas=0, sigmas_with_intra_block_xx=0, candidates=0, legal_sigmas_among_those=0, low_found=0)
found = []
for sig in itertools.permutations(range(6)):
    stats['sigmas'] += 1
    lab = rel(base, sig); g = [(idx(t), idx(s)) for t, s in lab]
    xx = [(t, s) for t, s in g if t < 6 and s < 6]
    assert len(xx) == 3
    cand = [j for j in range(6) if any(j in e and e[0] // 2 == e[1] // 2 for e in xx)]
    if not cand: continue
    stats['sigmas_with_intra_block_xx'] += 1; stats['candidates'] += len(cand)
    legal = search(g, ()) is not None
    stats['legal_sigmas_among_those'] += legal
    for j in cand:
        o = search(g, (j,)) if legal else None
        if o is not None:
            stats['low_found'] += 1; found.append(dict(sigma=sig, x=j, order=o))
    print(sig, 'cand', cand, 'legal', legal, 'found', len(found), '%.1fs' % (time.time() - t0), flush=True)
print(json.dumps(stats), 'seconds %.1f' % (time.time() - t0), flush=True)
print('FOUND', found[:5])
Path('relabel_fast.json').write_text(json.dumps(dict(stats=stats, found=found)))
