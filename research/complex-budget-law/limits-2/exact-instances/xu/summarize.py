"""expXU screen summary (stdlib, < 1 s): re-reads the logged screen results on the exact q = 9 instance (W = 166) and checks
the counts and nullities quoted in ../README.md.  It does not recompute the screens (they need python-flint, expXN/expQR
family builders and the upstream checker; see the README).
  GRS-F2 phase frames: screen_m{4..8}.json, 23 families x 3 char-2 draws, every nullity 0           -> 345 ranks full
  two-class theorem:   twoclass.json (35 partitions) + twoclass_cut.json (6 cut-aligned), no root of chi(M_S) in F_{2^d}, d <= 15
  rational rays:       ray_results.jsonl, 184 rays at m = 2..8, gcd degree 0 (all three draws and every pair), t = 1 nullity 0
Prints PASS_EXPXU_SUMMARY."""
import json, os
H = os.path.dirname(os.path.abspath(__file__))
ld = lambda f: json.load(open(os.path.join(H, f)))
ranks = 0
for m in range(4, 9):
    fams = {k: v for k, v in ld(f'screen_m{m}.json').items() if isinstance(v, dict) and 'nullity' in v}
    assert len(fams) == 23, (m, len(fams))
    for k, v in fams.items():
        assert v['nullity'] == [0, 0, 0] and v['Wm'] == 166 * m, (m, k, v)
        ranks += len(v['nullity'])
print(f'GRS-F2 screen: {ranks} ranks over 23 families x m = 4..8 x 3 draws, all full (nullity 0)')
tc = list(ld('twoclass.json').values()); cut = list(ld('twoclass_cut.json').values())
assert len(tc) == 35 and len(cut) == 6
assert all(r == [] for v in tc for r in v['roots_in_F2^d (d, multiplicity) per G draw'])
assert all(r == [] for v in cut for r in v['roots'])
print(f'two-class theorem: {len(tc) + len(cut)} partitions, no subfield root of chi(M_S) for d <= 15 on draw 1')
rays = [json.loads(l) for l in open(os.path.join(H, 'ray_results.jsonl'))]
assert len(rays) == 184 and {r['m'] for r in rays} == set(range(2, 9)) and all(r['W'] == 166 for r in rays)
assert all(r['gcd_deg_all'] == 0 and r['gcd_deg_pairs'] == [0, 0, 0] and r['nullity_t1_draw0'] == 0 for r in rays)
print(f'rational rays: {len(rays)} rays at m = 2..8, every gcd trivial, t = 1 nullity 0')
assert ranks == 345
print('PASS_EXPXU_SUMMARY')
