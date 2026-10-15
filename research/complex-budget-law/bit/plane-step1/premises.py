"""Verify the class-P premises for every pairwise-odd common-r set of a given family.
F1: the masked line target on L' = L \\ S forces X = sum y z^T to QF (81 F2 unknowns, full rank).
F2: for a in S, the only functional w with u_b.w = 0 for all b in L' with (b,a) a required zero is w = 0
    (equivalently those u_b span F2^9); same set serves rows by symmetry of the requirement.
usage: premises.py FAMILY  (FAMILY = all | 5_4 | s=K)"""
import core, json, sys, random, time
L = list(range(111)); U = core.U
def row(a, b): return sum((U[a][i] & U[b][j]) << (9 * i + j) for i in range(9) for j in range(9))
pairs = [(a, b) for a in L for b in L if a == b or core.G[a][b] % 2]
ROW = {p: row(*p) for p in pairs}
def F1(S):
    Sset = set(S); piv = {}; r = 0
    cand = [p for p in pairs if p[0] not in Sset and p[1] not in Sset]
    rnd = random.Random(len(cand)); rnd.shuffle(cand)
    for p in cand:
        v = ROW[p]
        while v:
            h = v.bit_length() - 1
            if h in piv: v ^= piv[h]
            else: piv[h] = v; r += 1; break
        if r == 81: return True
    return False
lab = [sum(U[t][i] << i for i in range(9)) for t in range(120)]
def F2(r, S):
    Sset = set(S)
    for a in S:
        vs = [lab[b] for b in L if b not in Sset and (core.G[a][b] % 2 or (core.G[a][b] == 0 and b != a and core.G[b][r] != 0))]
        if core.rank2(vs) != 9: return False
    return True
if __name__ == '__main__':
    fam = sys.argv[1]
    sets = json.load(open('sets_all.json'))
    if fam == '5_4': sets = [s for s in sets if len(s[1]) == 5 and s[2] == 4]
    elif fam.startswith('s='): sets = [s for s in sets if len(s[1]) == int(fam[2:])]
    t = time.time(); bad1 = bad2 = 0; seenS = {}
    for r, S, rk in sets:
        key = tuple(S)
        if key not in seenS: seenS[key] = F1(S)
        bad1 += not seenS[key]; bad2 += not F2(r, S)
    print(fam, 'sets', len(sets), 'distinct S', len(seenS), 'F1 failures', bad1, 'F2 failures', bad2, round(time.time() - t), 's')

