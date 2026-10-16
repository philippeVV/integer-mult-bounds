"""Sequence vs all-N (dense) exponent of the BGS amplification of the BH seed.
Uses BGS App. F recursion bound log m_i <= 4 c1^(2^i) * log(3 m_{i-1}) (edges), gap_i = (1+eps')^(2^i),
eps' from the alpha-extension: 1+eps' = (1+eps)/(1+alpha). Monotone gamma(N) => for N in [m_i, m_{i+1})
gamma(N) >= gap_i, so the all-N exponent is liminf log gap_i / log log m_{i+1}.
Works in log-log space to avoid overflow."""
import math
eps = 10/147; K = 314; k = 157; we = 1/785; hmax = 5
def run(alpha, label, cm=735+314):
    wl = alpha/(K*(1+eps)); c = (hmax*we + 2*wl)/min(wl, we); c1 = c + 1
    lg = math.log((1+eps)/(1+alpha))
    L = math.log(cm)                      # L_i = log m_i
    rows = []
    for i in range(1, 40):
        # log m_i <= 4 c1^(2^i) (log 3 + L_{i-1});  track log L
        logL = math.log(4) + (2**i)*math.log(c1) + math.log(math.log(3) + L) if L < 1e300 else None
        L = math.exp(logL) if logL < 690 else float('inf')
        rows.append((i, (2**i)*lg, logL))
        if L == float('inf'):
            # continue symbolically: logL_i ~ log4 + 2^i log c1 + logL_{i-1}
            prev = logL
            for j in range(i+1, 60):
                prev = math.log(4) + (2**j)*math.log(c1) + prev
                rows.append((j, (2**j)*lg, prev))
            break
    seq = rows[-2][1]/rows[-2][2]
    dense = rows[-3][1]/rows[-2][2]
    c2_printed = math.log(1+eps/2)/math.log(16*c1*c1) if label == 'BGS alpha' else lg/math.log(16*c1*c1)
    print(f"{label:12s} alpha={alpha:.5f} c1={c1:8.3f}  Thm1.1 c2={c2_printed:.5f}  recursion seq-limit={seq:.5f}  all-N(liminf)={dense:.5f}  (k=0.0017149: x{0.0017149/dense:.2f} of all-N)")
run((1+eps)/(1+eps/2)-1, 'BGS alpha')
best = max(((lambda a: a)(eps*i/2000), i) for i in range(1, 2000))
import io, contextlib
res = []
for i in range(1, 2000):
    a = eps*i/2000
    wl = a/(K*(1+eps)); c1 = (hmax*we + 2*wl)/min(wl, we) + 1
    res.append((math.log((1+eps)/(1+a))/(2*math.log(c1)), a))
r, a = max(res)
run(a, 'opt alpha')
