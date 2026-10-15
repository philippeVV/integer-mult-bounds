"""Independent re-check (no package imports): rank of any masked decoder N on a plane instance is > |S|+k.
Argument: S independent in Gamma_E => N[S,S] = I. If rank N = |S|+k then (Schur/Guttman)
  N = N[:,S] N[S,:] + U V^T,  U,V : n x k supported outside S.
For a not in S, s in S: N[a,s] != 0 only if E_a perp E_s (free pair).  Zero-pattern relaxation:
  diag a notin S : at least one nonzero product among {N[a,s]N[s,a]} u {U[a,i]V[a,i]}
  a != b notin S, (a,b) not free : number of nonzero products among {N[a,s]N[s,b]}, {U[a,i]V[b,i]} is != 1.
UNSAT => rank N >= |S|+k+1.
usage: rank_cert_indep.py TAR inst.json k s1,s2,...   (S given explicitly)"""
import sys, json, gzip, hashlib, tarfile, itertools
from pysat.solvers import Minisat22, Glucose4
TAR = sys.argv[1]; R = json.load(open(sys.argv[2])); R = (R[0] if isinstance(R, list) else R)['r']
k = int(sys.argv[3]); S = [int(x) for x in sys.argv[4].split(',')]
with tarfile.open(TAR) as t:
    raw = t.extractfile('candidate/strict/global/base/providers/gauge44/CANONICAL-WORD.json.gz').read()
assert hashlib.sha256(raw).hexdigest().startswith('16bc0ebb04b4be6c')
P = json.loads(gzip.decompress(raw))['ports']
P = [p if isinstance(p, int) else int(''.join(map(str, p)), 2) for p in P]
n = len(P)
dt = lambda a, b: bin(a & b).count('1') & 1
for q, r in zip(P, R): assert dt(q, q) == 1 and dt(r, q) == 0 and dt(r, r) == 1
free = [[a != b and all(dt(x, y) == 0 for x in (P[a], R[a]) for y in (P[b], R[b])) for b in range(n)] for a in range(n)]
assert len(set(S)) == len(S) and all(not free[a][b] for a in S for b in S), 'S must be pairwise non-free (independent)'
Sset = set(S); out = [a for a in range(n) if a not in Sset]
for Solver in (Minisat22, Glucose4):
    var = {}
    def V(*key):
        if key not in var: var[key] = len(var) + 1
        return var[key]
    cls = []
    def AND(x, y, key):  # t <-> x & y
        t_ = V(*key); cls.extend([[-t_, x], [-t_, y], [t_, -x, -y]]); return t_
    def products(a, b):
        ps = [AND(V('L', a, s), V('Rt', s, b), ('pr', a, s, b)) for s in S if free[a][s] and free[s][b]]
        ps += [AND(V('U', a, i), V('W', b, i), ('uv', a, b, i)) for i in range(k)]
        return ps
    for a in out:
        cls.append(products(a, a))
    for a, b in itertools.permutations(out, 2):
        if free[a][b]: continue
        ps = products(a, b)
        for p in ps: cls.append([-p] + [q for q in ps if q != p])
    with Solver(bootstrap_with=cls) as sol:
        res = sol.solve()
    print(Solver.__name__, 'k', k, '|S|', len(S), 'vars', len(var), 'clauses', len(cls), '=>', 'SAT (no conclusion)' if res else f'UNSAT: rank N >= {len(S) + k + 1}')
