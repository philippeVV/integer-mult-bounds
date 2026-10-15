"""5-subsets of the r-punctured 6-/7-sets whose own maxsav9 is 13: check maxsav9 <= 12 directly (n1 = 10 is covered
by monotonicity: no savings-19 configuration on the superset)."""
import json, itertools, sav9fast
P = list(range(111, 120))
sets = [s for s in json.load(open('sets_all.json')) if s[0] >= 111 and len(s[1]) >= 6]
bad = 0; n = 0; seen = set()
for r, S, rk in sets:
    if sav9fast.max_savings(P + S) <= 12: continue
    for sub in itertools.combinations(S, 5):
        if sub in seen: continue
        seen.add(sub); n += 1
        v = sav9fast.max_savings(P + list(sub))
        if v > 12: bad += 1; print('OPEN', r, sub, v)
print('checked 5-subsets', n, 'with maxsav9 > 12:', bad)
