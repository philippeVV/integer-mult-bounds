import sys, numpy as np, collections
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from e8 import *
R = np.array(ROOT, float) / 2
s = 0; t = next(b for b in range(NL) if IP[s][b] == 0 and b != s)
C2 = [b for b in range(NL) if b not in (s, t) and IP[b][s] != 0 and IP[b][t] != 0]
CL = {h: [a for a in range(NL) if dotq(LAB[a], h)] for h in range(1, 511)}
for h in (486, 504):
    O = CL[h]; out = [a for a in O if a not in C2]; inn = [a for a in O if a in C2]
    u = np.linalg.svd(R[out])[2][-1]
    print('h', h, '|O|', len(O), '|O minus C2|', len(out), '|O cap C2|', len(inn), 's,t in O', s in O, t in O,
          'u.r_s', round(float(u @ R[s]), 4), 'u.r_t', round(float(u @ R[t]), 4), '|u|', 1.0)
    print('   u vs (r_s +- r_t)/2:', round(float(abs(u @ (R[s] + R[t])) / np.linalg.norm(R[s] + R[t])), 4), round(float(abs(u @ (R[s] - R[t])) / np.linalg.norm(R[s] - R[t])), 4))
    print('   r_b.u on O cap C2 values', collections.Counter(round(float(R[b] @ u) * np.sqrt(2), 3) for b in inn), ' on C2 minus O', collections.Counter(round(float(R[b] @ u) * np.sqrt(2), 3) for b in C2 if b not in O))
print('labels: h486 ^ h504 =', 486 ^ 504, 'q_s', LAB[s], 'q_t', LAB[t], 'q_s^q_t', LAB[s] ^ LAB[t])
from maxdef import span2
rk = len(span2([LAB[b] for b in C2])).bit_length() - 1
std = [f for f in C2 if all(IP[f][b] != 0 for b in C2 if b != f)]
print('B = C2: rank q_B', rk, 'e9\' Def dim', 9 - rk, 'standard f (odd to all of B):', std)
print('486 in span(q_B)^perp?', all(dotq(486, LAB[b]) == 0 for b in C2), ' 504?', all(dotq(504, LAB[b]) == 0 for b in C2))
