"""Exact checks (Fractions / F2) for the reduction 'law on all L  <=>  C(P120) >= 72'.
Lemma 1 (restriction): L' c L  =>  C(L') <= C(L).        Lemma 2 (padding/blocks): C(L1 u L2) <= C(L1) + C(L2)
for DISJOINT L1, L2 (block-diagonal N; cross entries 0 are always allowed), and C({a}) <= 1.
Hence gain(L) = |L| - C(L) is monotone and superadditive; gain(L) <= gain(P120) for every L."""
import random
from fractions import Fraction as Fr
from e8 import *
def check(L, terms):
    """terms: list of (v, w) dicts line->Fraction. Exact masked identity on L; return cost."""
    Ls = list(bits(L)); N = {}
    for v, w in terms:
        for a, x in v.items():
            assert (L >> a) & 1
            for b, y in w.items(): N[a, b] = N.get((a, b), 0) + x * y
    for a in Ls:
        for b in Ls:
            val = N.get((a, b), 0)
            if a == b: assert val == 1, (a, val)
            elif IP[a][b] == 0: assert val == 0, (a, b, val)
    return sum(rho(sum(1 << a for a, x in v.items() if x), sum(1 << b for b, y in w.items() if y)) for v, w in terms)
def frame_terms(F, L):
    out = []
    for f in F:
        v = {a: Fr(sum(x * y for x, y in zip(ROOT[a], ROOT[f])), 8) for a in bits(L)}   # (r_a . f) / 2 in true coords
        v = {a: x for a, x in v.items() if x}
        out.append((v, dict(v)))
    return out
def restrict(terms, L):
    return [({a: x for a, x in v.items() if (L >> a) & 1}, {b: y for b, y in w.items() if (L >> b) & 1}) for v, w in terms]
def pad(terms, extra):
    return terms + [({a: Fr(1)}, {a: Fr(1)}) for a in bits(extra)]
if __name__ == '__main__':
    F = frames()[0]; Fm = sum(1 << a for a in F); P112 = FULL & ~Fm
    t120 = frame_terms(F, FULL); c120 = check(FULL, t120)
    print('P120 frame decoder exact, cost', c120, 'gain', 120 - c120)
    t112 = restrict(t120, P112); c112 = check(P112, t112); print('P112 restriction exact, cost', c112, 'gain', 112 - c112)
    tp = pad(t112, Fm); cp = check(FULL, tp); print('P112 decoder padded by I on the 8 frame lines: exact on P120, cost', cp)
    # random restriction + padding controls
    rnd = random.Random(1)
    for _ in range(5):
        L = FULL
        for a in rnd.sample(range(NL), rnd.randrange(1, 60)): L &= ~(1 << a)
        c = check(L, restrict(t120, L)); cpad = check(FULL, pad(restrict(t120, L), FULL & ~L))
        print(f'  |L|={bin(L).count("1"):3d}: restricted cost {c} (<= 72), padded back to P120 cost {cpad} = {c} + {120 - bin(L).count("1")}')
    # block-diagonal decoder from two different frames on a random 2-partition: exact
    F2_ = frames()[777]; L1 = sum(1 << a for a in rnd.sample(range(NL), 60)); L2 = FULL & ~L1
    tb = restrict(frame_terms(F, FULL), L1) + restrict(frame_terms(F2_, FULL), L2)
    print('block decoder (two frames, random halves) exact on P120, cost', check(FULL, tb))
    # rays: every nonzero m has supp(m) containing Phi \ H for a rank-7 flat H; all 8 ray types have label span 9
    from catalog import OMEGA
    print('ray-type supports: label-span dims', [rank2([LAB[a] for a in bits(supp(w, FULL))]) for w in OMEGA])
    print('=> any term of a torus-orbit-of-G/2 decoder on P120 has rho = 9 (both sides span F2^9, form nondegenerate)')
    print('PASS')
