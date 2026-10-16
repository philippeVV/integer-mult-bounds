"""BGS constant c1 = l/w + 1 for the Braverman-He seed (arXiv 2610.10108 Sec 4, App. A).
Rebuilds G (182 vertices, 735 edges), the 157 sessions B = R \\ (C u J), checks the
flow certificate (every session distance >= 5), takes the uniform dual w_e = 1/(157*5)
(objective 735/785 = 147/157), forms BGS's alpha-extension (App. C) and evaluates
c2 = ln(1+eps/2)/ln(16 c1^2) (BGS Thm 1.1 / App. F). Also scans alpha (same proof)."""
import math
from fractions import Fraction as Fr
from collections import deque
# F9 = F3[t]/(t^2+1); element r -> (r%3, r//3)
def mul(x, y):
    a, b = x; c, d = y
    return ((a*c - b*d) % 3, (a*d + b*c) % 3)
def add(x, y): return ((x[0]+y[0]) % 3, (x[1]+y[1]) % 3)
eta = [(r % 3, r//3) for r in range(9)]
Z, O = (0, 0), (1, 0)
T = [(Z, Z, O)] + [(Z, O, eta[r]) for r in range(9)] + [(O, eta[r], eta[s]) for r in range(9) for s in range(9)]
assert len(T) == 91
def dot(p, l):
    s = Z
    for a, b in zip(p, l): s = add(s, mul(a, b))
    return s
edges = []
for i in range(91):
    js = [j for j in range(91) if dot(T[i], T[j]) == Z]
    assert len(js) == 10
    for j in sorted(js): edges.append((i, 91+j))
assert len(edges) == 910
R = [7,8,9,30,31,38,39,56,57,60,68,75,79,84,88,94,102,103,117,122,127,128,141,151,154,155,157,163,
164,166,167,168,176,179,183,187,203,207,208,212,215,221,224,225,231,234,236,237,238,242,243,251,253,266,282,290,
300,309,318,321,324,325,330,333,338,341,349,380,386,387,389,396,398,400,401,406,407,409,411,412,418,425,426,428,
439,444,454,458,464,465,468,470,476,480,481,492,502,503,504,513,517,521,522,527,539,541,544,548,552,557,559,569,
579,583,585,588,592,593,598,602,603,606,607,612,619,620,632,634,643,646,670,672,673,674,682,685,688,694,695,701,
702,703,709,716,718,737,746,754,757,759,773,790,821,835,837,838,839,841,842,843,846,850,853,856,857,861,862,882,
884,885,887,900,901,904,907]
assert len(R) == 175 == len(set(R))
C = {102,163,164,341,386,406,476,513,862}; J = {103,349,470,517,603,606,674,716,861}
assert C <= set(R) and J <= set(R)
B = [q for q in R if q not in C | J]; assert len(B) == 157
Rs = set(R)
adj = [[] for _ in range(182)]
for q, (u, v) in enumerate(edges):
    if q not in Rs: adj[u].append(v); adj[v].append(u)
m = sum(len(a) for a in adj)//2; assert m == 735
def bfs(s):
    d = [-1]*182; d[s] = 0; dq = deque([s])
    while dq:
        x = dq.popleft()
        for y in adj[x]:
            if d[y] < 0: d[y] = d[x]+1; dq.append(y)
    return d
hs = [bfs(edges[q][0])[edges[q][1]] for q in B]
from collections import Counter
print("session hop distances:", sorted(Counter(hs).items()))
assert min(hs) >= 5
hmax, k = max(hs), len(B)
eps = Fr(157, 147) - 1
we = Fr(1, k*5)                    # uniform dual: sum_i dist_i = k*5*we = 1 (demand d = 1)
z = m*we; assert z == Fr(147, 157)
K = 2*k                            # leaves added: one per source and per sink (App. C)
cap = z*(1+eps)                    # leaf capacity z(D) d (1+eps) = 1
def c2_of(alpha):
    wl = alpha/(K*(1+eps))         # App. C leaf weight alpha/(k d (1+eps)), k = #leaves
    assert abs(K*cap*wl - alpha*z) < 1e-15    # leaves add exactly alpha*z(D)
    l = hmax*we + 2*wl             # largest session distance in D*
    w = min(wl, we)
    c1 = l/w + 1
    base = (1+eps)/(1+alpha)       # standard-form gap; amplified gap = base^(2^i)
    return float(c1), math.log(base)/math.log(16*c1*c1)
alpha_bgs = (1+eps)/(1+eps/2) - 1
print(f"eps = {eps} = {float(eps):.6f}; BGS alpha = {alpha_bgs} ; leaf weight / edge weight = {float(alpha_bgs/(K*(1+eps))/we):.5f}")
c1, c2 = c2_of(alpha_bgs)
print(f"BGS as stated: hmax={hmax}, c = l/w = {c1-1:.3f}, c1 = {c1:.3f}, c2 = {c2:.6f}, c2/0.0017149 = {c2/0.0017149:.3f}")
# alpha scan (same proof, alpha free in (0, eps))
best = max((c2_of(Fr(i, 10000)*eps)[1], i) for i in range(1, 10000))
a = Fr(best[1], 10000)*eps
print(f"alpha optimised: alpha/eps = {best[1]/10000:.4f}, c1 = {c2_of(a)[0]:.3f}, c2 = {best[0]:.6f}, /0.0017149 = {best[0]/0.0017149:.3f}")
# sensitivity: if a sharper (non-uniform) dual existed, c1 >= hmax-ish ; idealised lower bound on c1 with these leaves
for c1x in (2, 7, 20):
    print(f"  hypothetical c1={c1x}: c2 = {math.log(1+float(eps)/2)/math.log(16*c1x*c1x):.6f}")
