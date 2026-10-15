"""Extended finite premises: for classes O (h = 1..510) and lines f:
 m1 = min |O minus orth(f)| over O != E7A1(f);  m2 = min |O minus E7A1(f)| over O != E7A1(f);
 distribution of |O cap orth(s)|."""
import collections
from e8 import *
E7A1 = [ORTH[f] | (1 << f) for f in range(NL)]
CL = [sum(1 << a for a in range(NL) if dotq(LAB[a], h)) for h in range(1, 511)]
m1 = min(bin(O & ~ORTH[f]).count('1') for O in CL for f in range(NL) if O != E7A1[f])
m2 = min(bin(O & ~E7A1[f]).count('1') for O in CL for f in range(NL) if O != E7A1[f])
d = collections.Counter(bin(O & ORTH[s]).count('1') for O in CL for s in range(NL))
print('m1', m1, 'm2', m2, '|O cap orth(s)| distribution', sorted(d.items()))
