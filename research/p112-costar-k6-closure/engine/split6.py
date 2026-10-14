import pickle, itertools
d = pickle.load(open("certk6.pkl", "rb"))
filt = {tuple(A): (cs, dd) for A, cs, dd in d["filtered"]}
res = d["res"]; cov = d["cov"]; del d
comps = [c for c in itertools.product(range(1, 10), repeat=6) if sum(c) == 9]
jobs = []; opn = []
for A in cov:
    cs, dd = filt[tuple(A)]
    for c in comps:
        if c in dd: jobs.append((list(A), list(c), ("L1", dd[c]), []))
for x in res:
    if x[2] == "INFEASIBLE": jobs.append((list(x[0]), list(x[1]), x[3], x[4]))
    else: opn.append((x[0], x[1]))
print("jobs", len(jobs), "open", len(opn), "expected", len(cov) * len(comps))
n = 40
for i in range(n):
    pickle.dump(jobs[i::n], open(f"k6chunks/c{i:02d}.pkl", "wb"))
pickle.dump(opn, open("k6_open.pkl", "wb"))
