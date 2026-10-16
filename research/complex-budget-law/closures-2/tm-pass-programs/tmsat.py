"""Exact minimum step count of a k-tape binary oblivious TM realising a fixed
bit permutation (small transposes / address interchanges), XOR vs routing-only.

Model (TM-native, no frames):
  * k one-sided binary tapes, cells 0..L-1; tape 0 holds input bits x_0..x_{V-1}
    at cells 0..V-1, every other cell is 0 (blank = 0).  All heads start at cell 0.
  * finite control holds c data registers (bits), initially 0, plus an unbounded
    *program counter* (the program is non-uniform and oblivious: head moves and
    gate choices may depend on the step number, not on data).  This is a
    relaxation of a uniform TM, so its optimum is a LOWER bound for any real
    fixed machine on the same k, c.
  * one step: every head reads its cell; every tape writes one bit at its head and
    every register is updated, each as a function of the k read bits and c
    register bits; then every head moves -1/0/+1.
      mode 'xor'  : each written bit is any F2-linear combination of the k+c bits
      mode 'route': each written bit is a copy of ONE of the k+c bits, or 0
  * goal after T steps: tape `goal_tape` (or any tape if goal_tape is None)
    holds x_{perm[i]} at cell i for 0<=i<V.  Other cells/tapes/registers arbitrary.
Correctness is exact over all inputs: cells hold linear forms in F2^V.

The cost is T (= head moves of the busiest tape; every tape steps each time unit).
"""
import itertools, sys, time
from pysat.solvers import Solver
from pysat.card import CardEnc, EncType


class Enc:
    def __init__(self):
        self.n = 0
        self.cl = []
        self.T = self.new()
        self.cl.append([self.T])

    def new(self):
        self.n += 1
        return self.n

    def const(self, b):
        return self.T if b else -self.T

    def AND(self, a, b):
        if a == -self.T or b == -self.T: return -self.T
        if a == self.T: return b
        if b == self.T: return a
        y = self.new()
        self.cl += [[-y, a], [-y, b], [y, -a, -b]]
        return y

    def OR(self, xs):
        xs = [x for x in xs if x != -self.T]
        if any(x == self.T for x in xs): return self.T
        if not xs: return -self.T
        if len(xs) == 1: return xs[0]
        y = self.new()
        self.cl.append([-y] + xs)
        for x in xs: self.cl.append([y, -x])
        return y

    def XOR2(self, a, b):
        if a == -self.T: return b
        if b == -self.T: return a
        if a == self.T: return -b
        if b == self.T: return -a
        y = self.new()
        self.cl += [[-y, a, b], [-y, -a, -b], [y, -a, b], [y, a, -b]]
        return y

    def XOR(self, xs):
        r = -self.T
        for x in xs: r = self.XOR2(r, x)
        return r

    def MUX(self, s, a, b):  # s ? a : b
        if s == self.T: return a
        if s == -self.T: return b
        if a == b: return a
        y = self.new()
        self.cl += [[-s, -a, y], [-s, a, -y], [s, -b, y], [s, b, -y]]
        return y


def build(V, perm, k, c, T, L, mode, goal_tape=None, init=None):
    e = Enc()
    # content[j][p] = list of V literals (linear form coefficients)
    if init is None:
        init = [[[1 if (j == 0 and p < V and v == p) else 0 for v in range(V)]
                 for p in range(L)] for j in range(k)]
    cont = [[[e.const(init[j][p][v]) for v in range(V)] for p in range(L)] for j in range(k)]
    regs = [[-e.T] * V for _ in range(c)]
    pos = [[e.const(p == 0) for p in range(L)] for j in range(k)]
    nsrc = k + c
    for t in range(T):
        reads = []
        for j in range(k):
            reads.append([e.OR([e.AND(pos[j][p], cont[j][p][v]) for p in range(L)]) for v in range(V)])
        srcs = reads + regs
        outs = []
        for o in range(k + c):
            coef = [e.new() for _ in range(nsrc)]
            if mode == 'route':
                for cl in CardEnc.atmost(coef, 1, top_id=e.n, encoding=EncType.pairwise).clauses:
                    e.cl.append(cl)
                outs.append([e.OR([e.AND(coef[s], srcs[s][v]) for s in range(nsrc)]) for v in range(V)])
            else:
                outs.append([e.XOR([e.AND(coef[s], srcs[s][v]) for s in range(nsrc)]) for v in range(V)])
        # write
        newcont = []
        for j in range(k):
            newcont.append([[e.MUX(pos[j][p], outs[j][v], cont[j][p][v]) for v in range(V)] for p in range(L)])
        cont = newcont
        regs = outs[k:]
        # move
        newpos = []
        for j in range(k):
            mv = [e.new() for _ in range(3)]  # left, stay, right
            e.cl.append(mv)
            for a, b in itertools.combinations(mv, 2): e.cl.append([-a, -b])
            e.cl.append([-pos[j][0], -mv[0]])
            e.cl.append([-pos[j][L - 1], -mv[2]])
            row = []
            for p in range(L):
                terms = []
                if p - 1 >= 0: terms.append(e.AND(pos[j][p - 1], mv[2]))
                terms.append(e.AND(pos[j][p], mv[1]))
                if p + 1 < L: terms.append(e.AND(pos[j][p + 1], mv[0]))
                row.append(e.OR(terms))
            newpos.append(row)
        pos = newpos
    # goal
    tapes = range(k) if goal_tape is None else [goal_tape]
    sel = []
    for j in tapes:
        g = e.new(); sel.append(g)
        for i in range(V):
            for v in range(V):
                lit = cont[j][i][v]
                want = (v == perm[i])
                e.cl.append([-g, lit if want else -lit])
    e.cl.append(sel)
    return e


def solve(V, perm, k, c, T, L, mode, goal_tape=None, solver='cadical153', budget=None):
    e = build(V, perm, k, c, T, L, mode, goal_tape)
    with Solver(name=solver, bootstrap_with=e.cl) as s:
        r = s.solve()
    return r, e.n, len(e.cl)


def min_T(V, perm, k, c, L, mode, goal_tape=None, Tlo=1, Thi=200, log=True):
    """smallest T with a solution (monotone: an idle 'stay+keep' step is always available)."""
    T = Tlo
    while T <= Thi:
        t0 = time.time()
        r, nv, nc = solve(V, perm, k, c, T, L, mode, goal_tape)
        if log:
            print(f"   T={T:3d} {'SAT' if r else 'unsat'}  vars={nv} cls={nc} {time.time()-t0:.1f}s", flush=True)
        if r: return T
        T += 1
    return None


def rec_transpose(b, R=2):
    """R x R matrix of b-bit records, row-major -> transpose. Output cell i gets input perm[i]."""
    perm = []
    for i in range(R):
        for j in range(R):
            src = (j * R + i) * b
            perm += [src + z for z in range(b)]
    return perm


if __name__ == '__main__':
    pass
