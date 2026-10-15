"""Engine confirmation of Theorem 1 on every clique orbit (k = 2..8, non-A2): random B subset of C(S) with Clean(B) and random
admissible zero patterns (al = 0 / alp = 0 only on A cap B). Asserts max credit <= 9 (cost >= 72)."""
import json, random, sys
from relax import maxcredit
from orbits import clean
from e8 import *
orb = json.load(open('orbits.json')); rnd = random.Random(int(sys.argv[1])); N = int(sys.argv[2])
tot = 0; mx = 0
for k in map(str, range(2, 9)):
    for o in orb[k]:
        S = o['rep']; C = o['C']
        if o['rankS'] < len(S): print('k', k, 'A2 triple orbit skipped (not a kappa-family member, see README)'); continue
        n = 0; m = 0; tries = 0
        while n < N and tries < 20 * N:
            tries += 1
            p = rnd.random(); B = [c for c in C if rnd.random() < p]
            if not B: continue
            if clean(sum(1 << b for b in B)): continue      # non-empty list = not clean
            AB = [a for a in S if a in B]
            ZX = [a for a in AB if rnd.random() < 0.3]; ZY = [a for a in AB if rnd.random() < 0.3]
            v = maxcredit(S, B, ZX, ZY); m = max(m, v); n += 1
            assert v <= 9, (S, B, ZX, ZY, v)
        print('k', k, 'orbit', o['size'], '|C|', len(C), 'tested', n, 'max credit', m, flush=True); tot += n; mx = max(mx, m)
print('PASS' if mx <= 9 else 'FAIL', tot, 'cases, max credit', mx)
