"""Certificate (cheap facts) for the hub-relaxation obstructions; heavy exhaustive facts are re-derived by
exact56.py (maxports6=20, L6 with >=4 ports in L6^perp has <=12 ports), max5in7.py (5-dim inside a G26 edge
level <=10 ports, hyperplanes <=16) and betafree.py (no beta confines alpha_t to a proper subspace)."""
import json
from collections import Counter
P=json.load(open('instA.json'))['ports']; al=json.load(open('alpha.json')); n=112
nor=[1,7,25,31,121,127,97,103]; dot=lambda a,b:(a&b).bit_count()&1; idx=set(P)
S=lambda Q: frozenset(i for i,p in enumerate(P) if all(dot(p,q)==0 for q in Q))
ok=lambda c,m: print(('PASS ' if c else 'FAIL ')+m) or c
ok(all(dot(p,p)==1 for p in P) and len(idx)==112,'ports distinct, all odd (p.p=1; j=511 has j.p=1 for every port)')
ok(not (set(nor)&idx) and len(set(nor))==8,'8 distinct centre normals, none is a port => no U_q is a target hyperplane p_t^perp')
need={t:frozenset(int(s) for s in al[str(t)]) for t in range(n)}
ok(all(need[t]==S([P[t]]) for t in range(n)),'alpha_t support = all 52 ports in p_t^perp (spans it)')
m7=max(len(S([P[t],a])) for t in range(n) for a in range(1,512) if a!=P[t])
ok(m7==32,f'max ports in a 7-dim subspace of a target hyperplane = {m7} < 52  => every top needs >=2 contributing arrivals')
pair=Counter(len(S([P[t],P[u]])) for t in range(n) for u in range(n) if u!=t)
ok(max(pair)==26,f'ports in <p_t,p_u>^perp over all u!=t: {dict(pair)} (<=26)')
g=[(t,u) for t in range(n) for u in range(n) if u!=t and len(S([P[t],P[u]]))==26]
ok(all((P[t]^P[u]) not in idx for t,u in g),'26-port pair levels lie in exactly two tops (p_t+p_u never a port)')
print('2912 =',112*26,' > 224*12 =',224*12,' (dim-6 averaging contradiction for the helper-free design)')
