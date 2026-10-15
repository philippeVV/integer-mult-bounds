"""Class U closure for maximal sets with r punctured (where Odd_r labels have rank 8 and one extra rank per side is
possible).  Theta(s) = 64 + floor(8s/9) is the first failing cost.  n1 = 9: 72 - maxsav9 + 2 g(s) >= Theta.
n1 = 10: need no 10-term decoder with savings >= 81 - Theta + 2 g(s-1): s=5 -> 19 (ten19 exhaustive), s=6 -> 20 and
s=7 -> 21 impossible (savings 20 needs ten pair terms, but V n W must contain a non-pair term).  n1 >= 11 arithmetic."""
import json, sys, time
from collections import Counter
import sav9fast, ten19
P = list(range(111, 120))
def g(m): return m - 1 if m >= 2 else (1 if m == 1 else 0)
def theta(s): return 64 + (8 * s) // 9
sets = [s for s in json.load(open('sets_all.json')) if s[0] >= 111 and len(s[1]) >= 5]
t = time.time(); res = Counter(); mx = Counter()
for r, S, rk in sets:
    s = len(S); sav = sav9fast.max_savings(P + S)
    ok9 = 72 - sav + 2 * g(s) >= theta(s)
    need = 81 - theta(s) + 2 * g(s - 1)
    if need >= 20: ok10 = True
    else:
        assert need == 19
        ok10 = not ten19.exists19(P + S)[0]
    ok11 = all(6 * k + 2 * g(s - k + 9) >= theta(s) for k in range(11, 14))
    res[(s, ok9, ok10, ok11)] += 1; mx[(s, sav)] += 1
    if not (ok9 and ok10 and ok11): print('OPEN', r, S, sav, flush=True)
print('results', dict(res), 'sav9 distribution', dict(mx), round(time.time() - t), 's', flush=True)
