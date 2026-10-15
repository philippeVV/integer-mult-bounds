"""Independent replay verifier for co-star infeasibility certificates.

Does NOT import model.py / engine*.py.  Rebuilds E8 from scratch (doubled coordinates:
integer vectors y, all-even or all-odd, sum = 0 mod 4, |y|^2 = 8), orders the 120 lines
as in PR #1320, and replays each certificate with plain Python sets.

Lemma checks (all over any field of characteristic != 2, zeros allowed):
  L0: some live r has no group alive on both row and column at r.
  L1: orthogonal clique Q with sum_g min(c_g, |Q & R_g|, |Q & C_g|) < |Q|.
  L2: Q orthogonal clique, |Q| = c_g, all Q x Q entries owned by g alone; r orthogonal
      to all of Q with (r,q) [resp. (q,r)] owned by g alone for all q  -> kill row [col].
  L2rank: orthogonal clique of size c_g+1, all entries owned by g alone.
  L3: (c_g = 2) odd walk x..y in owned graph; r orth x,y with owned entries -> kill.
  L3odd: odd closed walk in the owned graph of a capacity-2 group.
Usage: verify_cert.py cert3.pkl cert4.pkl ...   (pickles produced by run_cert.py)
"""
import sys, pickle, itertools

# ---------- independent E8 ----------
def e8_doubled():
    roots = []
    for y in itertools.product(range(-2, 3), repeat=8):
        if sum(v * v for v in y) != 8:
            continue
        if len(set(v % 2 for v in y)) != 1:
            continue
        if sum(y) % 4 != 0:
            continue
        roots.append(y)
    return roots

ROOTS = e8_doubled()
assert len(ROOTS) == 240
def canon(y):
    # line representative: first nonzero coordinate positive
    for v in y:
        if v != 0:
            return y if v > 0 else tuple(-w for w in y)
LINES = sorted(set(canon(y) for y in ROOTS))
assert len(LINES) == 120
def order_key(v):
    if all(abs(x) == 1 for x in v):           # half-integral: after integral, product order
        return (1, tuple(0 if x > 0 else 1 for x in v[1:]))
    nz = [i for i, x in enumerate(v) if x]
    i, j = nz
    return (0, (i, j, 0 if v[j] > 0 else 1))
LINES.sort(key=order_key)
IP = lambda a, b: sum(x * y for x, y in zip(LINES[a], LINES[b])) // 4
DELETED = {0, 1, 26, 27, 44, 45, 54, 55}
assert [LINES[i] for i in sorted(DELETED)] == [
    (2, 2, 0, 0, 0, 0, 0, 0), (2, -2, 0, 0, 0, 0, 0, 0), (0, 0, 2, 2, 0, 0, 0, 0), (0, 0, 2, -2, 0, 0, 0, 0),
    (0, 0, 0, 0, 2, 2, 0, 0), (0, 0, 0, 0, 2, -2, 0, 0), (0, 0, 0, 0, 0, 0, 2, 2), (0, 0, 0, 0, 0, 0, 2, -2)]
assert LINES[56] == (1,) * 8 and LINES[87] == (1, 1, -1, -1, -1, -1, -1, -1)
import hashlib, json
# anchor the line order to PR #1320's published source_lines_sha256
assert hashlib.sha256(json.dumps([list(v) for v in LINES], sort_keys=True, separators=(',', ':')).encode()).hexdigest() \
    == "1c71b1394110e20a71a73e50a4b7b23e3ffce25e8bd3c27e04a94b738fff798e", "line order differs from PR1320"
LIVE = [i for i in range(120) if i not in DELETED]
ORTH = {(a, b) for a in range(120) for b in range(120) if a != b and IP(a, b) == 0}
def costar(alpha):
    return {r for r in LIVE if r == alpha or (alpha, r) in ORTH}
CS = {a: costar(a) for a in range(120)}

# ---------- replay ----------
class Fail(Exception):
    pass

def check(cond, msg):
    if not cond:
        raise Fail(msg)

def replay(normals, caps, final, steps):
    """positions in certificates are live positions; convert to source indices"""
    G = len(normals)
    check(len(set(normals)) == G and all(c >= 1 for c in caps) and sum(caps) == 9, "bad assignment")
    R = [set(CS[a]) for a in normals]
    C = [set(CS[a]) for a in normals]
    src = lambda p: LIVE[p]
    def alone(r, s, g):
        return all(not (r in R[h] and s in C[h]) for h in range(G) if h != g)
    def owned_vertex(x, g):
        return alone(x, x, g)
    def fzero(x, y):
        # M[x,y] forced zero: masked orthogonal, or no group alive at (row x, column y)
        return (x, y) in ORTH or all(not (x in R[h] and y in C[h]) for h in range(G))
    def owned_edge(x, y, g):
        return (x, y) in ORTH and alone(x, y, g) and alone(y, x, g)
    for st in steps:
        kind, g = st[0], st[1]
        side, r = st[2], src(st[3])
        if kind == "L2G":
            gam = list(g)
            Q = [src(q) for q in st[4]]
            Qs = set(Q)
            check(len(Qs) == len(Q) and all((p, q) in ORTH for p in Q for q in Q if p != q), "L2G clique")
            check(len(set(gam)) == len(gam) and sum(caps[h] for h in gam) == len(Q), "L2G tight")
            out = [h for h in range(G) if h not in gam]
            check(all(not (Qs & R[h] and Qs & C[h]) for h in out), "L2G block ownership")
            check(r not in Qs and all((r, q) in ORTH for q in Q), "L2G orth")
            if side == "R":
                check(all(not (r in R[h] and Qs & C[h]) for h in out), "L2G row ownership")
                for h in gam:
                    R[h].discard(r)
            else:
                check(all(not (r in C[h] and Qs & R[h]) for h in out), "L2G col ownership")
                for h in gam:
                    C[h].discard(r)
            continue
        if kind == "L3GZ":
            gam = list(g)
            check(sum(caps[h] for h in gam) == 2 and len(set(gam)) == len(gam), "L3GZ capacity")
            out = [h for h in range(G) if h not in gam]
            def ownG(a, b):
                return all(not (a in R[h] and b in C[h]) for h in out)
            x, y = src(st[4]), src(st[5])
            walk = [src(w) for w in st[6]]
            check(walk[0] == x and walk[-1] == y and (len(walk) - 1) % 2 == 1, "L3GZ walk parity")
            check(all(ownG(w, w) for w in walk), "L3GZ walk vertices")
            check(all(fzero(walk[i], walk[i + 1]) and fzero(walk[i + 1], walk[i]) and ownG(walk[i], walk[i + 1])
                      and ownG(walk[i + 1], walk[i]) and walk[i] != walk[i + 1] for i in range(len(walk) - 1)), "L3GZ walk edges")
            if side == "R":
                check(r != x and r != y and fzero(r, x) and fzero(r, y) and ownG(r, x) and ownG(r, y), "L3GZ row premise")
                for h in gam:
                    R[h].discard(r)
            else:
                check(r != x and r != y and fzero(x, r) and fzero(y, r) and ownG(x, r) and ownG(y, r), "L3GZ col premise")
                for h in gam:
                    C[h].discard(r)
            continue
        if kind == "L2GZ":
            gam = list(g)
            X = [src(q) for q in st[4]]
            Xs = set(X)
            check(len(Xs) == len(X) and all(fzero(p, q) for p in X for q in X if p != q), "L2GZ Z-clique")
            check(len(set(gam)) == len(gam) and sum(caps[h] for h in gam) == len(X), "L2GZ tight")
            out = [h for h in range(G) if h not in gam]
            check(all(not (Xs & R[h] and Xs & C[h]) for h in out), "L2GZ block ownership")
            check(r not in Xs, "L2GZ r in X")
            if side == "R":
                check(all(fzero(r, x) for x in X), "L2GZ row forced zero")
                check(all(not (r in R[h] and Xs & C[h]) for h in out), "L2GZ row ownership")
                for h in gam:
                    R[h].discard(r)
            else:
                check(all(fzero(x, r) for x in X), "L2GZ col forced zero")
                check(all(not (r in C[h] and Xs & R[h]) for h in out), "L2GZ col ownership")
                for h in gam:
                    C[h].discard(r)
            continue
        if kind == "L2":
            Q = [src(q) for q in st[4]]
            check(len(Q) == caps[g] and len(set(Q)) == len(Q), "L2 size")
            check(all(owned_vertex(q, g) for q in Q), "L2 diag ownership")
            check(all(owned_edge(p, q, g) for p in Q for q in Q if p != q), "L2 clique ownership")
            check(r not in Q and all((r, q) in ORTH for q in Q), "L2 orth")
            if side == "R":
                check(all(alone(r, q, g) for q in Q), "L2 row ownership")
                R[g].discard(r)
            else:
                check(all(alone(q, r, g) for q in Q), "L2 col ownership")
                C[g].discard(r)
        elif kind == "L3":
            check(caps[g] == 2, "L3 cap")
            x, y = src(st[4]), src(st[5])
            walk = [src(w) for w in st[6]]
            check(walk[0] == x and walk[-1] == y and (len(walk) - 1) % 2 == 1, "L3 walk parity")
            check(all(owned_vertex(w, g) for w in walk), "L3 walk vertices")
            check(all(owned_edge(walk[i], walk[i + 1], g) for i in range(len(walk) - 1)), "L3 walk edges")
            check((r, x) in ORTH and (r, y) in ORTH, "L3 orth")
            if side == "R":
                check(alone(r, x, g) and alone(r, y, g), "L3 row ownership")
                R[g].discard(r)
            else:
                check(alone(x, r, g) and alone(y, r, g), "L3 col ownership")
                C[g].discard(r)
        else:
            raise Fail("unknown step")
    kind = final[0]
    if kind == "L0":
        r = src(final[1])
        check(all(not (r in R[g] and r in C[g]) for g in range(G)), "L0 not violated")
    elif kind == "L1":
        Q = [src(q) for q in final[1]]
        check(len(set(Q)) == len(Q) and all((p, q) in ORTH for p in Q for q in Q if p != q), "L1 not a clique")
        tot = sum(min(caps[g], len(set(Q) & R[g]), len(set(Q) & C[g])) for g in range(G))
        check(tot < len(Q), "L1 not violated")
    elif kind == "L3GZodd":
        gam = list(final[1])
        check(sum(caps[h] for h in gam) == 2 and len(set(gam)) == len(gam), "L3GZodd capacity")
        out = [h for h in range(G) if h not in gam]
        def ownG(a, b):
            return all(not (a in R[h] and b in C[h]) for h in out)
        walk = [src(w) for w in final[2]]
        check(walk[0] == walk[-1] and (len(walk) - 1) % 2 == 1, "L3GZodd parity")
        check(all(ownG(w, w) for w in walk), "L3GZodd vertices")
        check(all(fzero(walk[i], walk[i + 1]) and fzero(walk[i + 1], walk[i]) and ownG(walk[i], walk[i + 1])
                  and ownG(walk[i + 1], walk[i]) and walk[i] != walk[i + 1] for i in range(len(walk) - 1)), "L3GZodd edges")
    elif kind == "L1Z":
        X = [src(q) for q in final[1]]
        check(len(set(X)) == len(X) and all(fzero(p, q) for p in X for q in X if p != q), "L1Z not a Z-clique")
        tot = sum(min(caps[g], len(set(X) & R[g]), len(set(X) & C[g])) for g in range(G))
        check(tot < len(X), "L1Z not violated")
    elif kind == "L2rank":
        g = final[1]
        Q = [src(q) for q in final[2]]
        check(len(Q) == caps[g] + 1, "L2rank size")
        check(all(owned_vertex(q, g) for q in Q) and all(owned_edge(p, q, g) for p in Q for q in Q if p != q), "L2rank own")
    elif kind == "L3odd":
        g = final[1]
        check(caps[g] == 2, "L3odd cap")
        walk = [src(w) for w in final[2]]
        check(walk[0] == walk[-1] and (len(walk) - 1) % 2 == 1, "L3odd parity")
        check(all(owned_vertex(w, g) for w in walk), "L3odd vertices")
        check(all(owned_edge(walk[i], walk[i + 1], g) for i in range(len(walk) - 1)), "L3odd edges")
    else:
        raise Fail("unknown final")
    return True

# ---------- coverage of the family by orbit reps ----------
def refl_perm(b):
    vb = LINES[b]
    key = {}
    for i, v in enumerate(LINES):
        key[v] = i
        key[tuple(-x for x in v)] = i
    out = []
    for v in LINES:
        ip = sum(x * y for x, y in zip(v, vb))
        out.append(key[tuple(x - (ip // 4) * y for x, y in zip(v, vb))])
    return out

def stabilizer_generators():
    Rf = [refl_perm(b) for b in range(120)]
    gens = []
    seen = set()
    for b in sorted(DELETED):          # single reflections in the deleted frame roots
        g = Rf[b]
        if {g[i] for i in DELETED} == DELETED and tuple(g) not in seen:
            seen.add(tuple(g)); gens.append(g)
    for b in range(120):               # double reflections preserving the frame
        for c in range(120):
            g = [Rf[b][Rf[c][i]] for i in range(120)]
            if {g[i] for i in DELETED} == DELETED and g != list(range(120)) and tuple(g) not in seen:
                seen.add(tuple(g)); gens.append(g)
    # sanity: each generator preserves |inner product| (it is induced by an isometry)
    for g in gens[:: max(1, len(gens) // 20)]:
        for a in range(0, 120, 11):
            for b in range(120):
                assert abs(IP(g[a], g[b])) == abs(IP(a, b))
    return gens

def covering_sets(k):
    full = set(LIVE)
    out = []
    for A in itertools.combinations(range(120), k):
        U = set()
        for a in A:
            U |= CS[a]
        if len(U) == 112:
            out.append(A)
    return out

if __name__ == "__main__":
    import time
    t0 = time.time()
    gens = stabilizer_generators()
    print("stabilizer generators (frame reflections + F0-preserving double reflections):", len(gens))
    total_cases = 0
    for path in sys.argv[1:]:
        res = pickle.load(open(path, "rb"))
        k = len(res[0][0])
        reps = sorted(set(x[0] for x in res))
        comps = {}
        for x in res:
            comps.setdefault(x[0], set()).add(tuple(x[1]))
        # every rep must come with all compositions of 9 into k positive parts
        allc = {c for c in itertools.product(range(1, 10), repeat=k) if sum(c) == 9}
        assert all(comps[r] == allc for r in reps), "missing compositions"
        ok = 0
        for x in res:
            normals, caps, status, final, steps = x[0], x[1], x[2], x[3], x[4]
            if status != "INFEASIBLE":
                print("NOT CLOSED:", normals, caps)
                continue
            replay(list(normals), list(caps), final, steps)
            ok += 1
        total_cases += ok
        print(path, ": k =", k, "reps", len(reps), "certificates verified", ok, "/", len(res))
        if k <= 4:
            fam = covering_sets(k) if k == 3 else None
            if k == 4:
                # covering 4-sets via pruning (same definition, independent code)
                fam = []
                full = set(LIVE)
                for a in range(120):
                    for b in range(a + 1, 120):
                        U2 = CS[a] | CS[b]
                        for c in range(b + 1, 120):
                            miss = full - U2 - CS[c]
                            if len(miss) > 60:
                                continue
                            for d in range(c + 1, 120):
                                if miss <= CS[d]:
                                    fam.append((a, b, c, d))
            idx = {s: i for i, s in enumerate(fam)}
            par = list(range(len(fam)))
            def find(z):
                while par[z] != z:
                    par[z] = par[par[z]]
                    z = par[z]
                return z
            for i, s in enumerate(fam):
                for g in gens:
                    j = idx[tuple(sorted(g[a] for a in s))]
                    ri, rj = find(i), find(j)
                    if ri != rj:
                        par[ri] = rj
            roots = {find(idx[r]) for r in reps}
            uncovered = sum(1 for i in range(len(fam)) if find(i) not in roots)
            print("   covering %d-sets: %d ; orbits %d ; reps given %d ; sets not in a rep orbit: %d"
                  % (k, len(fam), len({find(i) for i in range(len(fam))}), len(reps), uncovered))
    print("total verified certificates:", total_cases, " time %.1fs" % (time.time() - t0))
