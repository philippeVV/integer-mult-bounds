"""Class U (label supports on line ports unless one side of the term lies inside S; anchored sides arbitrary):
for r punctured the anchored-row (and -column) mechanism can add one extra rank via lambda_r (Odd_r labels have rank 8),
so the non-T1 bound weakens to 2*g(s-kappa), g(m) = m-1 (m>=2), 1 (m=1), 0 (m<=0).
Fail threshold Theta(s) = 64 + floor(8s/9).  Checks per maximal set with r in P:
  n1 = 9 : 72 - maxsav9(P u S) + 2 g(s) >= Theta
  n1 = 10: no 10-term label decoder with savings >= 81 - Theta + 2 g(s-1)  (exact F2 enumeration)
  n1 >= 11: 6 n1 + 2 g(s - n1 + 9) >= Theta (arithmetic)."""
import json, sys, time
from collections import Counter
import sav9fast, small_s, small12
P = list(range(111, 120))
def g(m): return m - 1 if m >= 2 else (1 if m == 1 else 0)
def theta(s): return 64 + (8 * s) // 9
def need10(s): return 81 - theta(s) + 2 * g(s - 1)
def check10(S, need):
    if need >= 21: return 'trivial'
    if need == 20: return small_s.check_set(S, 4)[0]           # ten pairs
    if need == 19: return small_s.check_set(S, 3)[0]           # nine pairs + one cost-7 (or ten pairs)
    if need == 18: return small12.check(S, 2)[0]
    if need == 17: return small12.check(S, 1)[0]
    raise ValueError(need)
if __name__ == '__main__':
    sizes = list(map(int, sys.argv[1].split(',')))
    sets = [s for s in json.load(open('sets_all.json')) if s[0] >= 111 and len(s[1]) in sizes]
    for s in sizes:
        for k in range(11, 14): assert 6 * k + 2 * g(s - k + 9) >= theta(s)
    t = time.time(); res = Counter(); worst = {}
    for r, S, rk in sets:
        s = len(S); sav = sav9fast.max_savings(P + S)
        ok9 = 72 - sav + 2 * g(s) >= theta(s)
        r10 = check10(S, need10(s))
        ok10 = r10 in ('none', 'nosol', 'trivial')
        res[(s, ok9, ok10)] += 1
        worst[s] = max(worst.get(s, 0), sav)
        if not (ok9 and ok10): print('OPEN', r, S, sav, r10, flush=True)
    print(dict(res), 'max sav9 by s', worst, 'need10', {s: need10(s) for s in sizes}, round(time.time() - t), 's', flush=True)
