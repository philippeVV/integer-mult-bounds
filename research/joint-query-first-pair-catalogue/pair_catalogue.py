#!/usr/bin/env python3
"""Exact test: does any alternative first route pair (F_a, B_{s,t}) of the QQ16
frame catalogue avoid the selected-observation loss?

Domain D: every replacement of the first pair F0/B_{1,2} by (F_a, B_{s,t}),
0<=a<=8, 1<=s<t<=8 (9*28 = 252 pairs), using the frozen logical maps (1),(3)
of CONSERVED-QQ16.md with routes D_u = rho_u(p), rho_0 = I, and the retained
ordinary tail / observation set of QQ16-SELECTED-WITNESS.md.

Transported witness: h_{s,t,k} = (48I - B_s* B_s) B_t* e_k.
Usage: pair_catalogue.py ARCHIVE WITNESS_JSON
"""
import gzip, json, sys, tarfile, time
from fractions import Fraction as Q
from hashlib import sha256
from io import BytesIO

A = [0, 7, 20, 27, 36, 46, 57, 111]
T = [21, 37, 47, 76, 81, 87, 92, 110]
U = [p for p in range(112) if p not in A + T]
MEMBERS = {
    "slots/local/CANONICAL-WORD.json.gz": "16bc0ebb04b4be6c0ade74badc6a0c0cdba4215eda86e4ff9f13abf8ecebd219",
    "slots/local/LOCAL-PROGRAMS.json.gz": "241fde27d4e43f5a9cb6b128f9f76c5131a53d09f0371676163f914d0b0fef3a",
    "t16/global/base/work/complex-nine-native/DATA-PHASES.json": "c5916ff312678f12aba403847643fbbc29f643acec60bdf836fa56cd3e8251fa",
}
ARCHIVE_SHA = "22c3f4dce4b897a1d423d1388801933e288b6a6b1906cb15e1b16a5eb68db754"


def need(c, m):
    if not c:
        raise SystemExit("FAIL: " + m)


def load(path):
    raw = open(path, "rb").read()
    need(sha256(raw).hexdigest() == ARCHIVE_SHA, "archive pin")
    out = {}
    with tarfile.open(fileobj=BytesIO(raw), mode="r:gz") as ar:
        for name, pin in MEMBERS.items():
            body = ar.extractfile(name).read()
            need(sha256(body).hexdigest() == pin, "member pin " + name)
            out[name.rsplit("/", 1)[1]] = json.loads(gzip.decompress(body) if name.endswith(".gz") else body)
    return out


def query_columns(local):
    word = local["CANONICAL-WORD.json.gz"]
    rows = [{i: Q(1)} if i < 112 else {} for i in range(712)]
    for ev in local["LOCAL-PROGRAMS.json.gz"]["forward"]:
        if ev == ["phase", "centers"]:
            break
        if ev[0] != "scalar":
            continue
        g = ev[1]
        reg = lambda t: t[1] + (704 if t[0] == "C" else 0)
        tg, op = reg(g[1]), g[0]
        if op == "clear":
            rows[tg] = {}
        elif op == "copy":
            rows[tg] = rows[reg(g[2])].copy()
        elif op == "scale":
            rows[tg] = {j: v * Q(g[2], g[3]) for j, v in rows[tg].items()}
        else:
            need(op == "add", "scalar alphabet")
            for j, v in rows[reg(g[2])].items():
                val = rows[tg].get(j, 0) + Q(g[3], g[4]) * v
                if val:
                    rows[tg][j] = val
                else:
                    rows[tg].pop(j, None)
    else:
        need(False, "centre cut")
    R = [[rows[c].get(p, Q(0)) for p in range(112)] for c in word["centers"]]
    need(all(v.denominator == 1 for r in R for v in r), "integer R")
    cols = [tuple(int(R[k][p]) for k in range(8)) for p in range(112)]
    gram = lambda cs: [[sum(v[i] * v[j] for v in cs) for j in range(8)] for i in range(8)]
    diag = lambda x: [[x * (i == j) for j in range(8)] for i in range(8)]
    need(all(sum(x * x for x in c) == 4 for c in cols), "norm four")
    need(gram([cols[p] for p in U]) == diag(48), "ordinary Gram 48I")
    need(gram([cols[p] for p in A]) == diag(4) and gram([cols[p] for p in T]) == diag(4), "sections")
    return cols


def act(cols, v):
    r = 0
    while v:
        lo = v & -v
        r ^= cols[lo.bit_length() - 1]
        v ^= lo
    return r


def mul(L, Rt):
    return tuple(act(L, c) for c in Rt)


def stage_route(basis, stage, coords):
    """Faithful restriction of rho_stage(p) to the invariant coordinate list coords
    (bits 0..8 plus external blocks of the stages involved), re-indexed compactly."""
    pos = {b: i for i, b in enumerate(coords)}
    comp = lambda v: sum(1 << pos[i] for i in range(73) if v >> i & 1)
    cols = []
    for b in coords:
        v = 1 << b
        if stage:
            for i, ax in enumerate(basis):
                sw = ax ^ (1 << (9 + 8 * (stage - 1) + i))
                if (v & sw).bit_count() % 2:
                    v ^= sw
        need(all(v >> i & 1 == 0 or i in pos for i in range(73)), "invariant restriction")
        cols.append(comp(v))
    return tuple(cols)


def coords_for(stages):
    c = list(range(9))
    for s in sorted(set(x for x in stages if x)):
        c += list(range(9 + 8 * (s - 1), 17 + 8 * (s - 1)))
    return c


def addm(tab, g, M, c):
    row = tab.setdefault(g, [0] * 64)
    for i in range(64):
        row[i] += c * M[i]


def outer(u, v):
    return [u[i] * v[j] for i in range(8) for j in range(8)]


def witness_and_queries(cols, bases, s, t, probes):
    """Literal group-ring h_{s,t} (all 8 k at once) and B_u h for u in probes.
    h is stored per port as {g: 8-vector over k}; B_u h as {g: 8x8}."""
    coords = coords_for([s, t] + list(probes))
    ID = tuple(1 << i for i in range(len(coords)))
    rho = {u: [stage_route(bases[p], u, coords) for p in range(112)] for u in set([s, t] + list(probes))}
    for u in rho:
        for p in U:
            need(mul(rho[u][p], rho[u][p]) == ID, "route involution")
    h = {}
    for p in U:
        row = {}
        # 48 * B_t* e_k : field p = r_p[k] g_t(p)^{-1}
        g = rho[t][p]
        row.setdefault(g, [0] * 8)
        for k in range(8):
            row[g][k] += 48 * cols[p][k]
        # - B_s* B_s B_t* e_k : field p = sum_q r_p.r_q r_q[k] g_s(p)^{-1} g_s(q) g_t(q)^{-1}
        for q in U:
            d = sum(a * b for a, b in zip(cols[p], cols[q]))
            if d:
                g = mul(rho[s][p], mul(rho[s][q], rho[t][q]))
                r = row.setdefault(g, [0] * 8)
                for k in range(8):
                    r[k] -= d * cols[q][k]
        h[p] = {g: v for g, v in row.items() if any(v)}
    out = {}
    for u in probes:
        tab = {}
        for p in U:
            for g, v in h[p].items():
                addm(tab, mul(rho[u][p], g), outer(cols[p], v), 1)
        out[u] = {g: m for g, m in tab.items() if any(m)}
    return h, out, ID


def tail(wjson):
    src = json.load(open(wjson))
    x, y, n = [0] * 4, [1, 0, 0, 0], 0
    for ev in src["trajectory_after_B1"]:
        if ev["kind"] == "bridge":
            d, s_, c = ev["dst"], ev["src"], ev["coefficient"]
            x[d] += c * x[s_]
            y[s_] -= c * y[d]
        elif ev["kind"] == "local":
            st = ev["stage"]
            if st % 2 == 0:
                need(not any(x), "later F source zero")
            else:
                need(y[st // 2] == 0, "later B source zero")
        else:
            need(ev["kind"] == "terminal", "known event")
            x, y = y, [-v for v in x]
        need([v for pr in zip(x, y) for v in pr] == ev["ordinary"], "trajectory agreement")
        n += 1
    return n, x, y


def main():
    t0 = time.time()
    local = load(sys.argv[1])
    cols = query_columns(local)
    ph = local["DATA-PHASES.json"]
    bases = [tuple(r["idle_children"][0]["orthogonal_axes"][:8]) for r in ph["recipes"]]
    word = local["CANONICAL-WORD.json.gz"]
    for p, r in enumerate(ph["recipes"]):
        loc = (r["q"],) + bases[p]
        need(r["port"] == p and loc[0] == word["ports"][p] and all(0 < v < 512 for v in loc), "chart")
        need(all((loc[i] & loc[j]).bit_count() % 2 == (i == j) for i in range(9) for j in range(9)), "orthonormal")
    res = {"pairs": {}, "selected_witness_h12": {}}
    I2112 = [2112 * (i == j) for i in range(8) for j in range(8)]
    # Part 1: transported witness for every (s,t); B_s h = 0 and B_t h, B_a h for all a.
    for s in range(1, 9):
        for t in range(s + 1, 9):
            # relative-route distinctness (full 73-bit faithful routes)
            full = coords_for([s, t])
            rel = [mul(stage_route(bases[p], t, full), stage_route(bases[p], s, full)) for p in U]
            distinct = len(set(rel)) == 96
            h, Bq, ID = witness_and_queries(cols, bases, s, t, [s, t])
            need(not Bq[s], f"B_{s} h_{s},{t} = 0")
            idc = Bq[t].get(ID)
            need(idc == I2112, f"identity coeff B_{t} h = 2112 I for ({s},{t})")
            res["pairs"][f"{s},{t}"] = {"relative_routes_distinct": distinct, "Bs_h_zero": True,
                                        "Bt_h_identity_coefficient": 2112, "Bt_h_group_terms": len(Bq[t]),
                                        "h_terms": sum(len(v) for v in h.values())}
            print(f"(s,t)=({s},{t}) distinct={distinct} Bs h=0, Bt h id=2112I, terms={len(Bq[t])}  {time.time()-t0:.0f}s", flush=True)
    # Part 2: the selected witness h_{1,2} under every other current chart s.
    for u in range(0, 9):
        h, Bq, ID = witness_and_queries(cols, bases, 1, 2, [u])
        res["selected_witness_h12"][str(u)] = {"Bu_h_zero": not Bq[u], "group_terms": len(Bq[u]),
                                               "identity_coefficient_diag": (Bq[u].get(ID) or [0] * 64)[0]}
        print(f"selected h_(1,2): B_{u} h zero={not Bq[u]} terms={len(Bq[u])}", flush=True)
    n, x, y = tail(sys.argv[2])
    res["tail"] = {"events": n, "final_x": x, "final_y": y}
    res["status"] = "PASS_QQ16_PAIR_CATALOGUE_WITNESS"
    json.dump(res, open(sys.argv[3], "w"), indent=1, sort_keys=True)
    print(res["status"], f"{time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
