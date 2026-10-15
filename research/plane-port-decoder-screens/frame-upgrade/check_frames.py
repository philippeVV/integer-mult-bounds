#!/usr/bin/env python3
r"""Independent checker (stdlib only): frame upgrade with a shared r on the SELECTED labels.
Ports: pinned canonical word in research/community-frontier1094-audit/selected-replay.tar.gz (sha256 16bc0ebb04b4be6c..).
E8 signs: plane-port-decoder-screens/lower/e8map.json, re-validated here (q.q' = <x,x'>+1 mod 2, Parseval through F0).
For every orthogonal 8-frame S of live ports and every r with r.q_a=0 (a in S), r.r=1, checks:
 (1) T := {b: r.q_b=1} \ S; every b in T is odd to exactly 4 frame ports and E_a not perp q_b for all of them;
     an 8-subset S'' of T is pairwise even, so U = S u S'' is independent in the new free graph: N[U,U]=I_16,
     hence rank N >= 16 and every decoder has >= 16 terms of positive cost.
 (2) rank_Q G[S,T] = rank_Q G[T,S] = 8.  Additive repair N = G/2 + X needs X[S,S]=0, X[S,T]=-G[S,T]/2, X[T,S]=-G[T,S]/2
     (all forced), so rank X >= 16 and the repair costs >= 16 (every positive-cost term costs >= 1).
 (3) construction I_S (+) (G/2)[S^c] (8 private plane terms of cost 2 + the 8 frame terms restricted to S^c) is a valid
     decoder of cost 80 (dC = 16), checked exactly in Fractions on a sample, cost on all.
"""
import gzip, hashlib, json, os, tarfile, itertools, sys, random
from fractions import Fraction as Q
TAR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../community-frontier1094-audit/selected-replay.tar.gz')
EMAP = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '../lower/e8map.json')
with tarfile.open(TAR) as t:
    raw = t.extractfile('candidate/strict/global/base/providers/gauge44/CANONICAL-WORD.json.gz').read()
assert hashlib.sha256(raw).hexdigest().startswith('16bc0ebb04b4be6c')
P = json.loads(gzip.decompress(raw))['ports']; P = [p if isinstance(p, int) else int(''.join(map(str, p)), 2) for p in P]
n = len(P); dot = lambda a, b: bin(a & b).count('1') & 1
em = {int(k): v for k, v in json.load(open(EMAP)).items()}
ip = lambda u, v: sum(x * y for x, y in zip(u, v)) // 4
X = [em[q] for q in P]; F0 = sorted(set(em) - set(P))
G = [[ip(X[a], X[b]) for b in range(n)] for a in range(n)]
assert all(G[a][a] == 2 and all(a == b or dot(P[a], P[b]) == (G[a][b] + 1) % 2 for b in range(n)) for a in range(n))
c = [[ip(X[a], em[f]) for f in F0] for a in range(n)]
assert all(sum(c[a][k] * c[b][k] for k in range(8)) == 2 * G[a][b] for a in range(n) for b in range(n))
D = [[a for a in range(n) if c[a][k]] for k in range(8)]
def red(vs):
    B = []
    for v in vs:
        for b in B: v = min(v, v ^ b)
        if v: B.append(v)
    return B
def prank(U, V):
    bu, bv = red(U), red(V); return len(red([sum(dot(x, y) << j for j, y in enumerate(bv)) for x in bu]))
def qrank(M):
    M = [[Q(x) for x in r] for r in M]; rk = 0
    for col in range(len(M[0])):
        p = next((i for i in range(rk, len(M)) if M[i][col]), None)
        if p is None: continue
        M[rk], M[p] = M[p], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][col]: f = M[i][col] / M[rk][col]; M[i] = [x - f * y for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk
assert sum(prank([P[a] for a in D[k]], [P[a] for a in D[k]]) for k in range(8)) == 64
EV = [[a != b and dot(P[a], P[b]) == 1 for b in range(n)] for a in range(n)]
frames = []
def ext(cur, cand):
    if len(cur) == 8: frames.append(cur); return
    for x in cand: ext(cur + [x], [y for y in cand if y > x and EV[x][y]])
for a in range(n): ext([a], [y for y in range(a + 1, n) if EV[a][y]])
assert len(frames) == 1260
def even_clique8(T):
    out = None
    def bk(cur, cand):
        nonlocal out
        if out: return
        if len(cur) == 8: out = cur; return
        for i, x in enumerate(cand): bk(cur + [x], [y for y in cand[i + 1:] if EV[x][y]])
    bk([], T); return out
cases = 0; c1cost = set(); rng = random.Random(0); exact_checked = 0
for S in frames:
    rs = [r for r in range(512) if all(dot(r, P[a]) == 0 for a in S) and dot(r, r) == 1]
    assert len(rs) == 16
    for k in range(8):
        pass
    base = sum(prank([P[a] for a in D[k] if a not in S], [P[a] for a in D[k] if a not in S]) for k in range(8))
    c1cost.add(base + 16)
    for r in rs:
        E = {a: (P[a], r) for a in S}
        T = [b for b in range(n) if b not in S and dot(r, P[b]) == 1]
        for b in T:
            odd = [a for a in S if G[a][b] % 2]
            assert len(odd) == 4 and all(any(dot(x, P[b]) for x in E[a]) for a in odd)
        U2 = even_clique8(T); assert U2 is not None
        U = list(S) + U2
        nonfree = lambda a, b: any(dot(x, y) for x in (E.get(a) or (P[a],)) for y in (E.get(b) or (P[b],)))
        assert all(nonfree(a, b) for a, b in itertools.combinations(U, 2))
        assert all(G[a][b] == 0 or not (G[a][b] % 2) for a in S for b in S if a != b)  # X[S,S]=0 forced
        assert qrank([[G[a][b] for b in T] for a in S]) == 8 and qrank([[G[b][a] for a in S] for b in T]) == 8
        cases += 1
        if rng.random() < 0.003:   # exact Fraction check of the dC=16 construction
            terms = [({a: Q(1)}, {a: Q(1)}) for a in S] + [({a: Q(c[a][k], 4) for a in D[k] if a not in S}, {a: Q(c[a][k]) for a in D[k] if a not in S}) for k in range(8)]
            N = {}
            for d, t in terms:
                for a, x in d.items():
                    for b, y in t.items(): N[a, b] = N.get((a, b), 0) + x * y
            assert all(N.get((a, a), 0) == 1 for a in range(n))
            assert all(N.get((a, b), 0) == 0 for a in range(n) for b in range(n) if a != b and nonfree(a, b))
            cost = sum(prank([w for a in d for w in (E.get(a) or (P[a],))], [w for b in t for w in (E.get(b) or (P[b],))]) for d, t in terms)
            assert cost == 80; exact_checked += 1
print('frames', len(frames), '(frame, r) cases', cases)
print('(1) every case: S u S\'\' independent of size 16 => rank N >= 16, >= 16 positive-cost terms')
print('(2) every case: X[S,T], X[T,S] forced with rank_Q 8 each => additive repair rank >= 16 => cost >= 16')
print('(3) construction I_S (+) (G/2)[S^c] costs', sorted(c1cost), '(dC = 16); exact Fraction checks passed on', exact_checked, 'random cases')
print('PASS_FRAME_UPGRADE_CHECK')
