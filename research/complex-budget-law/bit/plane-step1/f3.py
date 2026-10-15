"""F3: labels of Odd_r = {b in L \\ S : b.r != 0} span F2^9 (then Z_a, Y_a are a-independent linear functions of the
anchored corrections restricted to Odd_r, which every a in S must keep at zero). Monotone under shrinking S."""
import core, json, time
L = list(range(111))
lab = [sum(core.U[t][i] << i for i in range(9)) for t in range(120)]
sets = json.load(open('sets_all.json'))
t = time.time(); bad = 0; mn = 99
for r, S, rk in sets:
    Sset = set(S); vs = [lab[b] for b in L if b not in Sset and core.G[b][r] != 0]
    mn = min(mn, len(vs))
    bad += core.rank2(vs) != 9
print('F3 failures', bad, 'of', len(sets), 'min |Odd_r|', mn, round(time.time() - t), 's')
