import json, sys, collections
from fractions import Fraction as Q
d = json.load(open(sys.argv[1]))
exf = collections.defaultdict(list); ent = collections.Counter()
for ch, par, f in d['req']:
    exf[ch].append({int(k): Q(v) for k, v in f.items()})
    if not isinstance(par, list): ent[par] += 1
def rank(vecs):
    basis = []
    for v in vecs:
        v = dict(v)
        for p, b in basis:
            if v.get(p):
                c = v[p]
                for k, x in b.items():
                    y = v.get(k, 0) - c*x
                    if y: v[k] = y
                    else: v.pop(k, None)
        if v:
            p = min(v); c = v[p]; basis.append((p, {k: x/c for k, x in v.items()}))
    return len(basis)
dim = lambda V: (int(V).bit_count()).bit_length() - 1
B = collections.Counter(); L = 0
for V, fs in exf.items():
    b = len(fs) if dim(V) == 1 else len(fs) - rank(fs)
    B[dim(V)] += b; L += b
print('L', L, sorted(B.items()))
# waste view: L = arrivals + sum_V (entrants - rank(exit forms))
entf = collections.Counter()
for ch, par, f in d['req']:
    if not isinstance(par, list): entf[par] += 1
A = sum(1 for ch, par, f in d['req'] if isinstance(par, list))
Wd = collections.Counter(); Rd = collections.Counter(); Ed = collections.Counter(); Nd = collections.Counter()
for V in set(exf) | set(entf):
    r = 0 if dim(int(V)) == 1 else rank(exf.get(V, []))
    if dim(int(V)) == 1: continue
    Wd[dim(int(V))] += entf[V] - r; Ed[dim(int(V))] += entf[V]; Rd[dim(int(V))] += r; Nd[dim(int(V))] += 1
print('arrivals', A, 'waste', sum(Wd.values()), 'check', A + sum(Wd.values()))
for k in sorted(Wd): print(' dim', k, 'pools', Nd[k], 'entrants', Ed[k], 'outrank', Rd[k], 'waste', Wd[k])
