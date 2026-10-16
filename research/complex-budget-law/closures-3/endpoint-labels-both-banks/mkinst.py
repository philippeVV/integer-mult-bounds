"""expQW mkinst.py: write a tiny both-bank one-block instance for bfs.c.

Model (CM, both banks): one block F2^h, standard dot form. Frames = nondegenerate subspaces. Registers:
X_b (source <p_b>, sink q_b^perp), Y_b (source <q_b>, sink p_b^perp), H helpers (source 0, sink F).
Moves: a register steps to a nested frame (climb free, drop costs the dimension lost); a transvection r += c*s
over GF(P) is allowed when r and s share a frame (= a gate). Target values: X_b = -y_b, Y_b = x_b, helpers = z.
A finished register must sit inside its sink (climbing to the sink is free). Data drops are allowed.
usage: python3 mkinst.py h P H 'p1,q1;p2,q2;...' out.txt    (vectors as ints, bit i = e_{i+1})
"""
import sys, itertools

h, P, H = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
ports = [tuple(int(x) for x in pq.split(',')) for pq in sys.argv[4].split(';')]
dot = lambda u, v: bin(u & v).count('1') & 1


def span(vs):
    S = {0}
    for v in vs:
        S |= {s ^ v for s in S}
    return frozenset(S)


subs = set()
for k in range(h + 1):
    for vs in itertools.combinations(range(1, 1 << h), k):
        S = span(vs)
        if len(S) == 1 << k:
            subs.add(S)
def nondeg(S):
    return all(any(dot(u, w) for w in S) for u in S if u)
FR = sorted((S for S in subs if nondeg(S)), key=lambda S: (len(S), sorted(S)))
idx = {S: i for i, S in enumerate(FR)}
dim = [len(S).bit_length() - 1 for S in FR]
perp = lambda S: frozenset(v for v in range(1 << h) if all(dot(v, u) == 0 for u in S))
for p, q in ports:
    assert dot(p, p) == 1 and dot(q, q) == 1 and dot(p, q) == 0
V = len(ports); n = 2 * V + H
start, sink = [], []
for p, q in ports:
    start.append(idx[span([p])]); sink.append(idx[perp(span([q]))])     # X
for p, q in ports:
    start.append(idx[span([q])]); sink.append(idx[perp(span([p]))])     # Y
for _ in range(H):
    start.append(idx[frozenset({0})]); sink.append(idx[span(range(1 << h))])
# target rows: register i -> (column, coef); inputs ordered x_1..x_V, y_1..y_V, z_1..z_H
tgt = [(V + b, P - 1) for b in range(V)] + [(b, 1) for b in range(V)] + [(2 * V + j, 1) for j in range(H)]
with open(sys.argv[5], 'w') as f:
    f.write(f'{len(FR)} {n} {P} {2 * V}\n')
    for i, S in enumerate(FR):
        f.write(f'{dim[i]} ' + ''.join('1' if S <= T else '0' for T in FR) + '\n')   # row i: S_i subset S_j
    f.write(' '.join(map(str, start)) + '\n' + ' '.join(map(str, sink)) + '\n')
    f.write(' '.join(f'{c} {v}' for c, v in tgt) + '\n')
print(f'h={h} P={P} V={V} H={H}: {len(FR)} frames, n={n} registers, sum u = {2*V}')
