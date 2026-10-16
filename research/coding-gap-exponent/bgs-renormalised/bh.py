"""Rebuild the Braverman-He seed G (arXiv 2610.10108 Sec 4 / App. A): 182 vertices, 735 unit edges,
157 sessions (s_q = earlier endpoint in pi, t_q = later). Taken from expRL/bh_c1.py plus pi."""
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
PI = [93,23,16,143,107,54,88,38,157,40,171,124,123,118, 90,32,15,146,114,179,181,39,128,0,119,83,21,6,
132,160,17,79,22,42,41,174,5,167,26,96,8,61, 168,101,75,43,102,127,99,12,30,57,147,92,138,159,
113,31,172,7,116,78,28,126,19,120,89,117,9,14, 122,56,72,115,153,47,65,103,74,137,97,10,110,154,
135,95,4,125,86,139,34,162,80,87,129,27,134,140, 36,13,148,149,82,37,176,91,73,1,51,141,105,109,
35,44,130,66,11,98,121,25,76,169,49,53,150,2, 104,29,158,69,180,81,52,106,112,131,164,63,71,170,
155,24,64,156,62,48,46,145,163,45,77,175,33,142, 18,178,166,144,152,55,67,70,177,85,111,20,136,84,
100,58,151,68,50,133,165,94,60,173,3,59,161,108]
assert sorted(PI) == list(range(182))
POS = {v: i for i, v in enumerate(PI)}
CORE = [edges[q] for q in range(910) if q not in Rs]          # 735 physical edges
SESS = []                                                     # (q, s_q, t_q)
for q in B:
    u, v = edges[q]
    s, t = (u, v) if POS[u] < POS[v] else (v, u)
    SESS.append((q, s, t))
