"""Lifts whose support A (= S, rigid side) is NOT pairwise odd: A contains an orthogonal pair {s, s'}; then
B subset C(A) subset C({s, s'}) = lines non-orthogonal to both (s, s' themselves are excluded: each is orth to the other).
Clean is inherited by subsets, so Clean(C({s,s'})) for one orthogonal pair (W(E8) is transitive on orthogonal root pairs)
closes all such lifts by Theorem 1 (expSM). Also: no pairwise-odd 9-set exists (|S| <= 8)."""
import sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from orbits import clean, cliques
from e8 import *
from cliques_lib import NO
s = 0; t = next(b for b in range(NL) if IP[s][b] == 0 and b != s)
C2 = [b for b in range(NL) if b not in (s, t) and IP[b][s] != 0 and IP[b][t] != 0]
bad = clean(sum(1 << b for b in C2))
print('orth pair', (s, t), '|C(A)|', len(C2), 'Clean bad classes', bad)
# transitivity check on orthogonal pairs: count orbits via generators (pairs as 2-sets)
from orbits import GEN
pairs = [(a, b) for a in range(NL) for b in range(a + 1, NL) if IP[a][b] == 0]
idx = {p: i for i, p in enumerate(pairs)}; par = list(range(len(pairs)))
def f(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for i, (a, b) in enumerate(pairs):
    for g in GEN:
        q = tuple(sorted((g[a], g[b]))); j = idx[q]; x, y = f(i), f(j)
        if x != y: par[x] = y
print('W-orbits of orthogonal line pairs:', len({f(i) for i in range(len(pairs))}))
print('pairwise-odd 9-sets:', len(cliques(9)))
# non-Clean B inside C({s,t}) forces B = C2 (both bad classes contain all of C2, nc2.py). Can a larger A keep C(A) = C2?
ext = [a for a in range(NL) if a not in (s, t) and all(IP[a][b] != 0 for b in C2 if b != a)]
print('lines a (not s,t) non-orth to all of C2:', ext, '-> A = {s,t} is the only non-clique support with non-Clean B'
      if not ext else '')
import numpy as np
from maxdef import span2
print('B = C2: rank q_B', len(span2([LAB[b] for b in C2])).bit_length() - 1, '; standard h\'_f:', [f for f in C2 if all(IP[f][b] != 0 for b in C2 if b != f)])
