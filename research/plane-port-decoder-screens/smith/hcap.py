"""EXACT caps H(rho,d) = max r_d of the 0/1 pairing matrix odd(P) x odd(Q) (odd vectors != j) over maximal pairs (P,Q),
P ⊇ Q^perp, rank(P,Q) = dim P + dim Q - 9 = rho, for rho = 9, 8 (all three hyperplane types, both orientations)."""
import sys
from smithnp import r_j
from smith import g
def dot(a, b): return bin(a & b).count('1') & 1
J = 511
ODD = [x for x in range(1, 512) if bin(x).count('1') & 1 and x != J]
def M(Pset, Qset): return [[dot(x, y) for y in Qset] for x in Pset]
def inperp(y): return [x for x in ODD if dot(x, y) == 0]
reps = {'even y': 3, 'odd y!=j': 1, 'y=j': J}   # 3 = e0+e1 even; 1 = e0 odd
D = range(0, 4)
out = {}
full = M(ODD, ODD)
out[(9, 'V,V')] = [r_j(full, d) for d in D]
for name, y in reps.items():
    Q = inperp(y); mm = M(ODD, Q)
    out[(8, 'V,' + name)] = [r_j(mm, d) for d in D]
for k, v in out.items(): print('rho', k[0], k[1], ' r_d for d=0..3:', v, ' g_d:', [g(d, k[0]) for d in D], flush=True)
