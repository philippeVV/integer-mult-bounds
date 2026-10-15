"""Top-level test (hub relaxation, no helpers): top p_t^perp has exactly 2 arrivals (t, s_t) with distinct L7's.
An L7 lying in only one top carries <= maxports(6-dim) = 20; a 'paired' L7 = <p_t,p_u>^perp carries <= its ports.
Need union = all 52 ports of p_t^perp. Also general: best union over ANY two L7 (port count only)."""
import json
from collections import Counter
P=json.load(open('instA.json'))['ports']; n=112; dot=lambda a,b:(a&b).bit_count()&1
S=lambda Q: frozenset(i for i,p in enumerate(P) if all(dot(p,q)==0 for q in Q))
ok=Counter(); anyb=Counter()
for t in range(n):
    full=S([P[t]]); assert len(full)==52
    pl=[S([P[t],P[u]]) for u in range(n) if u!=t]
    good=any(len(a|b)==52 for i,a in enumerate(pl) for b in pl[i+1:])
    ok[good]+=1
    hs=[S([P[t],a]) for a in range(1,512) if a!=P[t]]
    anyb[max(len(a|b) for i,a in enumerate(hs) for b in hs[i+1:])]+=1
print('targets with a covering pair of paired L7s:',dict(ok))
print('max |union| over ANY two L7s (port bound only):',dict(anyb))
