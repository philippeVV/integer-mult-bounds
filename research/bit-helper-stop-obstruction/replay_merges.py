"""For each of the 205 rank-1 stops build the natural merge word (stop gate (F1,h,[s]) moved to h's
next frame F2, inserted immediately before h's first F2 touch; C rebuilt as the literal H4 reverse),
then run: value-only 661-column endpoint, the semantic D/K boundary of the PR1311 bridge, the pinned
xcheck_p full replay (frames on), and the bridge's literal COPY/ERASE local bill."""
import gzip, json, sys, copy, importlib.util
from collections import Counter, defaultdict
ROOT = sys.argv[1]
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
xc = load('xc', ROOT+'/research/e8-bit-p111-tuned/code/xcheck_p.py')
br = load('br', ROOT+'/research/r439-abc-semantic-bridge/verify.py')
W = json.loads(gzip.decompress(open(ROOT+'/research/r439-abc-semantic-bridge/inputs/candidate.json.gz','rb').read()))
O = xc.Word(W); v=W['v']; n=2*v+W['R']
is_z = lambda r: not (v <= r < 2*v)
def rebuild_C(w):
    zg = [[O.FULL, g[1], list(g[2])] for sec in ('A','B') for g in w[sec] if is_z(g[1])]
    w['C'] = list(reversed(zg))
assert (lambda w: (rebuild_C(w), w['C'] == W['C'])[1])(copy.deepcopy(W)), 'C is literal reverse'
seq = [('D',i) for i in range(len(W['D']))]+[('A',i) for i in range(len(W['A']))]+[('B',i) for i in range(len(W['B']))]
touch = defaultdict(list)
for t,(ph,i) in enumerate(seq):
    f,d,src = W[ph][i]
    for r in [d]+src: touch[r].append((t,f,ph,i))
def owner_for(w):
    o = copy.copy(O); o.w = w; return o
# baseline sanity
base = owner_for(W); base.full(); bill0 = br.local_bill(W, base); N0 = br.boundary(W)
print('baseline', bill0['ADDs'], bill0['profile'], bill0['rank_mass'])
out = Counter(); first = {}
for h in range(2*v, n):
    p=[W['start'][h]]; ft={}
    for t,f,ph,i in touch[h]:
        if f!=p[-1]: p.append(f); ft[f]=(t,ph,i)
    if not (len(p)>=3 and O.fdim[p[1]]==1 and O.fdim[p[2]]==2): continue
    _,ph0,i0 = ft[p[1]]; _,ph2,i2 = ft[p[2]]
    assert ph0=='A' and ph2 in 'AB'
    w = copy.deepcopy(W)
    g = w['A'][i0]; s = g[2][0]
    new = [p[2], h, [s]]
    if ph2=='A':
        w['A'].insert(i2, new); del w['A'][i0]
    else:
        del w['A'][i0]; w['B'].insert(i2, new)
    rebuild_C(w)
    o = owner_for(w)
    # value-only endpoint over all 661 columns
    val,_,_ = o.replay(frames_on=False)
    vals_ok = all(val[r]==1<<r for r in range(n) if is_z(r)) and all(val[v+t]==(1<<(v+t))|(1<<t) for t in range(v))
    try: sem = br.boundary(w) == N0
    except ValueError as e: sem = 'reject:'+str(e)
    try: o.full(); xr = 'pass'
    except xc.Reject as e: xr = 'reject:'+str(e)
    except Exception as e: xr = 'reject:'+str(e)[:40]
    try: br.local_bill(w, o); lb='pass'
    except ValueError as e: lb='reject:'+str(e)
    key=(vals_ok, sem, xr, lb); out[key]+=1
    first.setdefault(key, (h, s, w['frames'][p[2]]))
for k,c in out.items(): print(c, k, 'e.g.', first[k])
