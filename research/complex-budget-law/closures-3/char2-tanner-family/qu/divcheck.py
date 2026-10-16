"""Check: random block orders (Singer-subgroup orbits, contiguous) -> exact full mcf |R| divisible by |H|."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import random, numpy as np
from baer import singer_D, G
n, D = singer_D(4); g = G(16, n, D, 59); rng = random.Random(5)
ok = True
for m in (21, 39, 91, 13):
    H = n // m; out = []
    for t in range(4):
        blocks = list(range(2 * m)); rng.shuffle(blocks)
        if t == 0: blocks = [x for k in range(m) for x in (k, m + k)]
        order = [v for b in blocks for v in ([b + m * s for s in range(H)] if b < m else [n + (b - m) + m * s for s in range(H)])]
        prof, pos = g.profile(order); sl, R = g.mcf(pos); out.append((R, R % H, sl))
    ok &= all(R % H == 0 for R, _, _ in out)
    print("m", m, "|H|", H, "(|R|, |R| mod |H|, slack):", out, flush=True)
print('PASS_EXPQU_DIV' if ok else 'FAIL_EXPQU_DIV')
