"""Write a BlockInvert-format promotion recipe (data only) for E1 -> E2 promotions with minimal-completion gate raises.
usage: make_e1_recipe.py WORD.json.gz cones.json OUT-RECIPE.json H:x [H:x ...]   (expB)
The recipe is consumed unchanged by research/gauge44-leaver8/construction/promo/build_promotion.py."""
import sys, json, gzip, hashlib, importlib.util
from collections import Counter
sys.path.insert(0, __import__('os').path.dirname(__file__))
import e1lib as L
s = importlib.util.spec_from_file_location('build_promotion', str(__import__('pathlib').Path(__file__).resolve().parents[3]) + '/research/gauge44-leaver8/construction/promo/build_promotion.py')
B = importlib.util.module_from_spec(s); s.loader.exec_module(B)
word, cones, out = sys.argv[1:4]
raw = gzip.decompress(open(word, 'rb').read()); w = json.loads(raw)
W = L.Word(w); W.sha = hashlib.sha256(raw).hexdigest(); W.prepare(); W.compute_cones(cones); W.ceilings()
sel = {}
for p in sys.argv[4:]:
    h, x = map(int, p.split(':')); sel[h] = L.add(W.start[h], L.mask_of_basis((x,)))
raised, dprof, legal, bad, touched = W.complete(sel)
assert legal, bad
gf = []
for gi in sorted(raised):
    sec, i, rs, U = W.ev[gi]
    old = w[sec][i][-1] if sec in ('A', 'B') else w['early_cleanup'][i][1]
    gf.append([sec, i, list(old), list(L.basis_of_mask(raised[gi]))])
r = dict(input_content_sha256=W.sha, entrances={str(h): list(L.basis_of_mask(m)) for h, m in sorted(sel.items())}, gate_frames=gf)
nw = B.apply(w, r)
data = json.dumps(nw, separators=(',', ':')).encode()
r['expected'] = dict(output_content_sha256=hashlib.sha256(data).hexdigest(), local_profile=nw['report']['local_profile'],
                     width_inventory=nw['report']['width_inventory'])
old = Counter({int(k): c for k, c in w['report']['local_profile'].items()}); new = Counter({int(k): c for k, c in nw['report']['local_profile'].items()})
d = Counter(new); d.subtract(old); d = Counter({k: c for k, c in d.items() if c})
assert d == dprof, (d, dprof)
json.dump(r, open(out, 'w'), indent=1)
print('recipe', out, 'entrances', len(sel), 'gate frames', len(gf), 'dprofile', dict(sorted(d.items())), 'profile', r['expected']['local_profile'],
      'widths', r['expected']['width_inventory'], 'sha', r['expected']['output_content_sha256'])
