"""rank of the natural decoder matrix N = G/2 - A_lost/2 = I + A_E/2 (A_E = E8 odd part kept on Gamma_E)."""
import json, sys
from gfamily import R, ip, rank_mod, n
from gamma import edges
for f in sys.argv[1:]:
    data = json.load(open(f)); runs = data if isinstance(data, list) else [data]
    for run in runs:
        r = run['r']; E = edges(r)
        N = [[2 if i == j else 0 for j in range(n)] for i in range(n)]  # 2N to stay integral
        for i, j in E: N[i][j] = ip(R[i], R[j]); N[j][i] = ip(R[j], R[i])
        print(f, 'edges', len(E), 'rank(I + A_E/2) =', rank_mod(N), flush=True)
