"""expQE exact single-column replacement bound for the x^4 code on schedule784's delay network.
Replace any one column of G_L (or of H_P) by an ARBITRARY v in F4096^8 (all rescalings, point moves, non-MDS codes):
dim W(v) = base_d + #{eta : v in V_eta} <= base_d + 273.  Checks base_d for all 32 columns, that the current code gives 772,
and that max_d base_d + 273 <= 781 < 784, so no single-column change makes the schedule exact.  (The RANSAC search of s2.py,
which found nothing above 772, is not rerun.)"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from lib import *
import colrep
from circ import gen_twisted
CP, CL = gen_twisted(Dp, 2)
GL = mul(CL, np.array([apow(9 * d) for d in Dp])[None, :]); HP = nullspace(CP)
bases = []; ok = True
for side in ['L', 'P']:
    for d in range(16):
        H, G = (HP, GL) if side == 'L' else (GL, HP)
        base, eqs = colrep.systems(H, G, d)
        cur = base + colrep.count(eqs, G[:, d])
        bases.append(base); ok &= cur == 772
        print(side, d, 'base', base, 'current dim W', cur, 'upper bound', base + 273, flush=True)
print('max base', max(bases), '-> any single-column replacement gives dim W <=', max(bases) + 273)
print('PASS_EXPQE_COLREP' if ok and max(bases) + 273 <= 781 else 'FAIL_EXPQE_COLREP')
