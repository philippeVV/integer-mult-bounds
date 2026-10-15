"""Stop-test probe of the |S| = 2 non-Clean locus: B = F u S with kap'_b = lam * (r_b . r_s1) on F (so the extra Y-line E1
exists), optionally kap'_s1 / al'_s1 tuned so that h'_s1 = E1 (merged, Def 2-dim). Exact member, exact deficient lines,
exact-rho gauge bound mincost (necessary conditions; uses the actual pairings of this point)."""
import json, random, sys, numpy as np
from kfam import RT, dot, check, deficient, Fr
from maxdef import mincost
from e8 import *
orb = json.load(open('orbits.json')); S = orb['2'][0]['rep']; C = orb['2'][0]['C']; F = [c for c in C if c not in S]
s1, s2 = S; rnd = random.Random(int(sys.argv[1]))
def rq(): return Fr(rnd.choice([-1, 1]) * rnd.randint(1, 9), rnd.randint(1, 9))
for merged in (False, True):
  for trial in range(int(sys.argv[2])):
    lam = rq(); X = []; Y = []
    for a in range(NL):
        ka = rq() if a in S else Fr(0)
        if a in F: kp = lam * dot(RT[a], RT[s1]); al = Fr(1); alp = (1 - ka * kp) / 2
        elif a in S:
            al = rq(); kp = rq()
            if merged and a == s1:   # E1 = (r_s1, -1/lam') vanishes on F with (kp_b = lam r_b.r_s1, alp_b = 1/2): r_b.r_s1/2 + kp_b g = 0 -> g = -1/(2 lam)
                g = -1 / (2 * lam)    # require 2 alp_s1 + kp_s1 g = 0 and 2 al alp + ka kp = 1
                # alp = -kp g / 2 ; 2 al (-kp g/2) + ka kp = 1 -> kp (ka - al g) = 1
                if ka == al * g: ka += 1
                kp = 1 / (ka - al * g); alp = -kp * g / 2
            else: alp = (1 - ka * kp) / (2 * al)
        else: kp = Fr(0); al = Fr(1); alp = Fr(1, 2)
        X.append([al * c for c in RT[a]] + [ka]); Y.append([alp * c for c in RT[a]] + [kp])
    assert check(X, Y)
    dX, kx = deficient(X); dY, ky = deficient(Y)
    def L(d): return [(np.array([float(x) for x in l]) / np.linalg.norm([float(x) for x in l]), h) for l, h in d.items() if l[0] != 'dim']
    mc = mincost(L(dX), L(dY)) if kx <= 1 and ky <= 1 else None
    print('merged' if merged else 'E1 only', 'X lines', sorted(len(h) for h in dX.values()), 'Y lines', sorted(len(h) for h in dY.values()),
          'kerdim', kx, ky, 'mincost', mc[0] if mc else None, flush=True)
