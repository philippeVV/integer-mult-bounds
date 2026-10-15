"""Monotonicity cover for non-maximal r-punctured sets of size 5 and 6 inside maximal 6- and 7-sets:
maxsav9 and the savings-19 configurations are monotone in P', so a maximal set S with maxsav9(P u S) <= 12 and no
savings-19 configuration certifies every 5-subset (class U, n1 = 9 and 10); maxsav9 <= 13 certifies every 6-subset."""
import json, time
from collections import Counter
import sav9fast, ten19
P = list(range(111, 120))
sets = [s for s in json.load(open('sets_all.json')) if s[0] >= 111 and len(s[1]) >= 6]
t = time.time(); res = Counter()
for r, S, rk in sets:
    sav = sav9fast.max_savings(P + S)
    no19 = not ten19.exists19(P + S)[0]
    res[(len(S), 'sav<=12' if sav <= 12 else ('sav<=13' if sav <= 13 else 'sav>13'), 'no19' if no19 else 'has19')] += 1
print(dict(res), round(time.time() - t), 's', flush=True)
