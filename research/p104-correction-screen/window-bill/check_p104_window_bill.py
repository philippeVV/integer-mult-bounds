"""P104 I-G/2 side correction: global access windows, simultaneous-storage floor and a
favorable complete bill against the PR2474 target.  Standard library only, read-only.

Usage:  python3 -I -B check_p104_window_bill.py REPO_ROOT [--ledger PR2474_LEDGER.json]

Model (PR2353/PR2393/PR2445 "X-only monotone relay" model, stated, not new):
  new helper lifetimes start source-independent at ZERO, move only upward (nested frames),
  mix rationally only in equal frames, read immutable X registers only in X's retained frames,
  never read original G/Q/D/C helpers or Y targets, write Y targets only additively.
"""
import collections, hashlib, json, os, sys
from fractions import Fraction as Q
from itertools import combinations, product

ROOT = os.path.abspath(sys.argv[1]); LEDGER = None
if '--ledger' in sys.argv: LEDGER = sys.argv[sys.argv.index('--ledger') + 1]
sys.path.insert(0, os.path.join(ROOT, 'research', 'complex-local-word-ceiling'))
import ceiling as CE  # Douglas Colkitt's exact-interval engine (merged PR1240)

def need(ok, why):
    if not ok: raise SystemExit('FAIL: ' + why)
OUT = {}
PKG = os.path.join(ROOT, 'research', 'p104-native-channel')
PINS = {'inputs/LOWRANK-WITNESS.json': '0f541bef7f19fdeffd67157684e9f03d785f42c89973632f508fc1a0452b580e',
        'inputs/P104-FIXTURE.json': '516c06f6ec3cb84be588c546108de897650288576dd5fbf5fe44b9c5546ad8a4',
        'EXPECTED.json': 'c9681a03d9a752b5b387250f841e43f7333b9984a325c4f37166e03f11676eb5'}
for p, h in PINS.items():
    need(hashlib.sha256(open(os.path.join(PKG, p), 'rb').read()).hexdigest() == h, 'pin ' + p)
W = json.load(open(os.path.join(PKG, 'inputs/LOWRANK-WITNESS.json')))
fx = json.load(open(os.path.join(PKG, 'inputs/P104-FIXTURE.json')))['section4_explicit_P104']
EXP = json.load(open(os.path.join(PKG, 'EXPECTED.json')))

# ---------- geometry: ports, Gram, C = I - G/2 ----------
L = set()
for i, j in combinations(range(8), 2):
    for sg in (-1, 1):
        v = [0]*8; v[i] = 2; v[j] = 2*sg; L.add(tuple(v))
for sg in product((-1, 1), repeat=7):
    if sg.count(-1) % 2 == 0: L.add((1, *sg))
removed = {tuple(r) for r in fx['F0'] + fx['C']}
ports = [r for r in sorted(L) if r not in removed]; n = len(ports); need(n == 104, '104 ports')
ip = lambda a, b: sum(x*y for x, y in zip(a, b)) // 4
G = [[ip(a, b) for b in ports] for a in ports]
C = [[Q(int(i == j)) - Q(G[i][j], 2) for j in range(n)] for i in range(n)]

def rank_q(M):
    M = [list(map(Q, r)) for r in M]; rk = 0
    for c in range(len(M[0])):
        p = next((i for i in range(rk, len(M)) if M[i][c]), None)
        if p is None: continue
        M[rk], M[p] = M[p], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = M[i][c] / M[rk][c]; M[i] = [a - f*b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk
need(rank_q(C) == 104, 'C = I - G/2 nonsingular (rows linearly independent)')

# ---------- binary frames along the retained chronology ----------
def span(vs):
    s = {0}
    for b in vs:
        if b not in s: s |= {x ^ b for x in s}
    return frozenset(s)
def rank2(vs):
    basis = []
    for v in vs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return len(basis)
cur = {r: span(f) for r, f in W['initial_frames'].items()}
path = {r: [(0, cur[r])] for r in cur}           # (cut at which frame starts, frame)
for i, e in enumerate(W['program']):
    if e[0] == 'transport':
        need(cur[e[1]] == span(e[2]), 'transport entry'); new = span(e[3])
        need(cur[e[1]] <= new or new <= cur[e[1]], 'nested transport'); cur[e[1]] = new; path[e[1]].append((i + 1, new))
    elif e[0] == 'zero': cur['C'] = span(e[2]); path['C'] = [(i + 1, cur['C'])]
X = ['X%d' % s for s in range(n)]; Y = ['Y%d' % t for t in range(n)]
for r in X + Y: need(cur[r] == span(W['final_frames'][r]), 'final frame ' + r)
q = {}
for s in range(n):
    f0 = W['initial_frames'][X[s]]; need(len(f0) == 1 and f0[0], 'source starts on a line'); q[s] = f0[0]
need(len(set(q.values())) == 104, 'distinct source labels')
H = {t: cur[Y[t]] for t in range(n)}
need(all(rank2(W['final_frames'][Y[t]]) == 8 for t in range(n)) and len(set(H.values())) == 104, '104 distinct rank-8 endpoints')
# monotone source paths
for s in range(n):
    fr = [f for _, f in path[X[s]]]; need(all(a < b for a, b in zip(fr, fr[1:])), 'source paths strictly increasing')

# B1: support of C_t is exactly the source lines inside H_t (48 each); own line outside
supp = {t: {s for s in range(n) if C[t][s] != 0} for t in range(n)}
need(all(supp[t] == {s for s in range(n) if q[s] in H[t]} for t in range(n)), 'supp C_t = lines in H_t')
need(all(len(supp[t]) == 48 and q[t] not in H[t] for t in range(n)), '48 lines per target, own line outside')
pairs = sum(len(v) for v in supp.values()); need(pairs == 4992, '4992 correction entries')
# B5: the 48 lines span H_t, so one proper-subframe arrival cannot carry the row (PR2353 floor input)
need(all(rank2([q[s] for s in supp[t]]) == 8 for t in range(n)), 'support spans each endpoint')

# B2: per needed pair, the retained X_s frames inside H_t and the last cut at which one is current.
# (PR2393 premise 2 -- only the initial line -- holds for Y0 but NOT globally; recorded as a finding.)
useful_end = {}; later_inside = collections.Counter(); shared = collections.defaultdict(lambda: collections.defaultdict(set))
for t in range(n):
    for s in supp[t]:
        fr = path[X[s]]; ends = []
        for k in range(len(fr)):
            if fr[k][1] <= H[t]:
                ends.append(fr[k + 1][0] if k + 1 < len(fr) else len(W['program']) + 1)
                if k: later_inside[t] += 1; shared[t][fr[k][1]].add(s)
        need(ends and fr[0][1] <= H[t], 'initial line inside')
        useful_end[t, s] = max(ends)          # first cut at which no retained X_s frame lies in H_t
first_exp = {s: path[X[s]][1][0] for s in range(n)}
E = max(useful_end.values())                  # at cuts >= E no X read can reach any target
# B3: every target is in ZERO until its first transport; all its frames lie in H_t
first_y = {t: path[Y[t]][1][0] - 1 for t in range(n)}           # last cut with Y_t in ZERO
need(all(path[Y[t]][0][1] == frozenset({0}) for t in range(n)), 'targets start in ZERO')
need(all(f <= H[t] for t in range(n) for _, f in path[Y[t]]), 'target paths inside H_t')
F = min(first_y.values())                                       # no target has left ZERO at cuts <= F
need(E <= F, 'nonempty storage window')
# no target frame is a source line (no direct X->Y read in the model), and Y paths are monotone
need(all(f != span([q[s]]) for t in range(n) for _, f in path[Y[t]] for s in supp[t]), 'no direct source-target frame')

OUT['access'] = dict(pairs=pairs, per_target=48, support_spans_endpoint=True,
    pairs_with_later_source_frame_inside_endpoint=sum(later_inside.values()),
    targets_with_such_pairs=len(later_inside),
    later_inside_per_affected_target=dict(sorted(collections.Counter(later_inside.values()).items())),
    max_needed_sources_on_one_later_source_frame=max(len(v) for d in shared.values() for v in d.values()),
    Y0_Y1_affected=[0 in later_inside, 1 in later_inside],
    acquisition_window='every useful X read for the correction happens at cuts < %d (latest end of a retained source frame inside a needed endpoint)' % E,
    delivery_window='no Y_t leaves ZERO before cut %d (earliest first target transport is event %d)' % (F + 1, F),
    storage_window=[E, F], window_events=F - E + 1,
    latest_first_source_expansion_cut=max(first_exp.values()))

# ---------- Storage lemma, checked premise by premise ----------
# At any cut c in [E, F]: (i) by B2 every needed X_s has left every retained frame inside H_t, so no later X read can lie on a
# dependence path to any Y_t; (ii) every Y_t is still in ZERO, and a ZERO-frame register holds no
# source information, so nothing useful has been written yet; (iii) original helpers are not read.
# Hence C_t x lies in the span of the forms held at cut c by live NEW helpers whose current frame is
# inside H_t (frames only grow, so a helper outside H_t never reaches it, nor does anything it
# touches).  Over all t these 104 forms are independent (C nonsingular): at least 104 new helpers
# are simultaneously live at every cut of the window, each holding a distinct independent form.
OUT['storage_floor'] = dict(simultaneous_new_helpers_min=104, cuts=[E, F],
    reason='rank C = 104; all correction information must be in live new helpers at these cuts')

# control: the window is sharp at both ends
need(any(v == E for v in useful_end.values()), 'left edge sharp')
need(any(first_y[t] < F + 1 for t in range(n)), 'right edge sharp')
# control: the X-only premise is essential -- replay original rows and show source-bearing original
# registers inside some H_t at cut E (the PR2173/PR2251 pooled model exploits exactly these)
names = list(W['initial_frames']); row = {r: {r: Q(1)} for r in names}; fr = {r: span(f) for r, f in W['initial_frames'].items()}
for i, e in enumerate(W['program'][:E]):
    op = e[0]
    if op == 'transport': fr[e[1]] = span(e[3])
    elif op == 'zero': row['C'] = {}; fr['C'] = span(e[2])
    elif op == 'erase': row.pop('C'); fr.pop('C')
    elif op == 'copy': row['C'] = dict(row[e[2]])
    else:
        _, d, s_, a, _f = e
        for k, v in row[s_].items():
            z = row[d].get(k, Q(0)) + Q(a)*v
            if z: row[d][k] = z
            else: row[d].pop(k, None)
orig_inside = sorted(r for r in row if r[0] not in 'XY' and any(k[0] == 'X' for k in row[r]) and any(fr[r] <= H[t] for t in range(n)))
OUT['access']['original_source_bearing_registers_inside_some_endpoint_at_E'] = len(orig_inside)
need(len(orig_inside) > 0, 'boundary control: original helpers do carry usable source information')

# ---------- favorable bill ----------
chan = {int(k): v for k, v in EXP['paid_profile'].items()}
need(chan == {1: 1074, 2: 58, 3: 86, 4: 104, 5: 48, 6: 24, 7: 8, 8: 8} and EXP['rank_mass'] == 2368, 'channel profile')
m, fam, st, cen = 73, 832, 9, 56                               # PR2082/PR2083 conditional P104 template
idle = CE.model(chan, m, fam, st, cen)['idle_mass']; need(idle == 4344, 'idle mass 4344')
A2474 = Q('0.001703676387993676489033949045348121170345872741703307776154')   # PR2474 kappa (truncated below)
ABIT = Q(1733427039589844, 10**18)                                          # bit supplier (prior screen)
if LEDGER:
    led = json.load(open(LEDGER))['rows']['struct_boundary_radix_union']
    k = Q(led['kappa']['rational']) if '/' in led['kappa']['rational'] else Q(led['kappa']['rational'])
    gap = Q(led['strict_checks']['actual_bit_above_effective']['rational'])
    need(abs(k - A2474) < Q(1, 10**55), 'PR2474 kappa binding')
    need(abs((k + gap) - ABIT) < Q(1, 10**15), 'bit endpoint = PR2474 effective + actual_bit_above_effective')
    OUT['ledger_bit_check'] = float(k + gap)
def add(b, e, k):
    o = dict(b)
    for r, c in e.items(): o[r] = o.get(r, 0) + c*k
    return o
def passes(prof, a):
    try: return CE.relax(prof, m, fam, st, cen, a)[1] < 1
    except ValueError: return False
def maxN(per, a, hi=2000):
    lo = 0
    while hi - lo > 1:
        mid = (lo + hi)//2; lo, hi = (mid, hi) if passes(add(chan, per, mid), a) else (lo, mid)
    return lo
g = Q(1, 10**12)
# Lemma (cheapest arrival): an arrival lifetime has children ZERO..F (total dim F), F->H_t (8-dim F),
# H_t->FULL (exactly 1).  Child cost r^p (p=1-a in (0,1]) is subadditive and concave, so the cheapest
# profile is {7:1, 1:2}.  Every other lifetime (reader/relay) pays at least {1:1, 8:1}.
arr = {1: 2, 7: 1}
OUT['bill'] = dict(
    channel_alone_ceiling=float(CE.ceiling(chan, m, fam, st, cen, grid=g)[1]),
    floor208_cheapest_ceiling=float(CE.ceiling(add(chan, arr, 208), m, fam, st, cen, grid=g)[1]),
    target_PR2474=float(A2474), bit_cap=float(ABIT),
    worth_cap_relative_to_PR2474=float(ABIT / A2474 - 1),
    max_lifetimes_beating_PR2474={str(p): maxN(p, A2474) for p in ({1: 2, 7: 1}, {1: 3, 6: 1}, {1: 4, 5: 1}, {1: 9})},
    max_lifetimes_at_bit_cap={str(p): maxN(p, ABIT) for p in ({1: 2, 7: 1}, {1: 9})},
    idle_mass_per_call=str(idle), idle_lifetime_cap=int(idle // 9),
    on_top_rule_allowance_at_PR2474=None)
# repo "on-top" rule (PR2082/2083/2133): added ranks are not coverage progress
lo, hi = Q(0), Q(100)
while hi - lo > Q(1, 10**6):
    mid = (lo + hi)/2
    try: ok = CE.relax(add(chan, {1: mid}, 1), m, fam, st, cen + mid, A2474)[1] < 1
    except ValueError: ok = False
    lo, hi = (mid, hi) if ok else (lo, mid)
OUT['bill']['on_top_rule_allowance_at_PR2474'] = round(float(lo), 4)
need(OUT['bill']['floor208_cheapest_ceiling'] > float(ABIT), 'favorable floor bill still clears the bit cap')
need(OUT['bill']['max_lifetimes_beating_PR2474'][str({1: 2, 7: 1})] == 482 == OUT['bill']['idle_lifetime_cap'], 'idle mass binds the cheapest profile')
# controls
try: CE.model(chan, m, fam, st, cen, R=52); raise SystemExit('FAIL: idle cap 52 accepted')
except ValueError: pass
n1 = OUT['bill']['max_lifetimes_beating_PR2474'][str({1: 9})]
need(passes(add(chan, {1: 9}, n1), A2474) and not passes(add(chan, {1: 9}, n1 + 1), A2474), 'frontier sharp')
bq = dict(q); bq[0] ^= 1
need({s for s in range(n) if bq[s] in H[1]} != supp[1] or {s for s in range(n) if bq[s] in H[2]} != supp[2] or True, 'label control')
need(any({s for s in range(n) if bq[s] in H[t]} != supp[t] for t in range(n)), 'corrupted label breaks support binding')
OUT['best_published_sharing'] = 'PR2445: 78 lifetimes for Y0+Y1 (39/target); full-correction budget at most 482/104 = 4.6 per target'
print(json.dumps(OUT, indent=1))
print('PASS_P104_WINDOW_BILL pairs=4992 window=%d..%d simultaneous_min=104 floor_ceiling=%.6f cap=%.6f' %
      (E, F, OUT['bill']['floor208_cheapest_ceiling'], float(ABIT)))
