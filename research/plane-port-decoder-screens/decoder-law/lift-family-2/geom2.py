import sys, numpy as np, collections
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from e8 import *
from maxdef import span2
R = np.array(ROOT, float) / 2
s = 0
E = [a for a in range(NL) if IP[a][s] == 0 and a != s]; N = [a for a in range(NL) if IP[a][s] != 0 and a != s]
sg = {b: np.sign(R[b] @ R[s]) for b in N}
p = {b: sg[b] * R[b] - R[s] / 2 for b in N}            # r_b (sign-fixed, r_b.r_s = 1) = r_s/2 + p_b
CL = {h: [a for a in range(NL) if dotq(LAB[a], h)] for h in range(1, 511)}
h0 = [h for h, O in CL.items() if not set(O) & set(E)]; print('h0', h0, sorted(CL[h0[0]]) == sorted(N))
hb = {}
for h, O in CL.items():
    OE = [a for a in O if a in E]
    if len(OE) == 36:
        # direction
        u = np.linalg.svd(np.vstack([R[OE], R[s]]))[2][-1]
        b = [b for b in N if abs(abs(p[b] @ u) / np.linalg.norm(p[b]) - 1) < 1e-9]
        hb[h] = b
        ON = [c for c in O if c in N]
        dist = collections.Counter(round(float(p[c] @ p[b[0]]), 3) for c in ON)
        print('h', h, 'dir weights', b, '|O cap N|', len(ON), 's in O', s in O, 'p_c.p_b distribution on O cap N', sorted(dist.items()))
        break
print('span of h0 + all h_b:', len(span2(h0 + [h for h, O in CL.items() if len([a for a in O if a in E]) == 36])).bit_length() - 1)
print('q_s', LAB[s], 'h0.q_s', dotq(h0[0], LAB[s]))
print('full p.p distribution', collections.Counter(round(float(p[c] @ p[N[0]]), 3) for c in N))
