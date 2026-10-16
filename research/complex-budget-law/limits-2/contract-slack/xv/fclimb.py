"""Hill-climb the order to minimise F(sigma) (min back arcs of an Eulerian orientation), starting from orient_NAME.json order.
usage: fclimb.py NAME steps seed"""
import sys, json, os, time, numpy as np
sys.argv, args = [sys.argv[0], sys.argv[1], '0', '1'], sys.argv
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'orient.py')).read().split('rng = np.random.default_rng(0)')[0])
steps, seed = int(args[2]), int(args[3]); rng = np.random.default_rng(seed)
d0 = json.load(open(f'orient_{name}.json')); order = np.array(d0['order'])
def Fo(order):
    key = np.empty(n); key[order] = np.arange(n); return F_of(key)
cur = Fo(order); t0 = time.time()
for s in range(steps):
    o2 = order.copy(); i = rng.integers(n); j = int(np.clip(i + rng.integers(-20, 21), 0, n - 1))
    x = o2[i]; o2 = np.delete(o2, i); o2 = np.insert(o2, j, x)
    r = Fo(o2)
    if (r[1], r[0]) <= (cur[1], cur[0]): order, cur = o2, r
    if s % 200 == 0: print(s, cur[0], cur[1], round(time.time() - t0), flush=True)
c, Fv, tail, head, back, pos = cur
if Fv < d0['F']:
    d0.update(cost=c, F=Fv, order=[int(x) for x in order], tail=tail.tolist(), head=head.tolist(), R=back.tolist())
    json.dump(d0, open(f'orient_{name}.json', 'w'))
print('final', name, c, Fv)
