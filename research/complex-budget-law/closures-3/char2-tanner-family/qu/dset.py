"""Planar Singer difference set of PG(2,2^e) via trace F_{2^{3e}} -> F_{2^e}; cached to dset_e.npy."""
import numpy as np, os
POLY = {2: (1 << 6) | 0b11, 4: (1 << 12) | (1 << 6) | (1 << 4) | 2 | 1, 6: (1 << 18) | (1 << 7) | 1}
def dset(e):
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"dset_{e}.npy")
    if os.path.exists(fn): return [int(x) for x in np.load(fn)]
    deg = 3 * e; poly = POLY[e]; q = 2 ** e; n = q * q + q + 1; ORD = 2 ** deg - 1
    EXP = [0] * ORD; LOG = [0] * (ORD + 1); v = 1
    for k in range(ORD):
        EXP[k] = v; LOG[v] = k; v <<= 1
        if v >> deg & 1: v ^= poly
    assert v == 1 and len(set(EXP)) == ORD
    def pw(y, k): return 0 if y == 0 else EXP[(LOG[y] * k) % ORD]
    D = [i for i in range(n) if (EXP[i] ^ pw(EXP[i], q) ^ pw(EXP[i], q * q)) == 0]
    assert len(D) == q + 1
    s = set(); [s.add((a - b) % n) for a in D for b in D if a != b]; assert len(s) == n - 1
    np.save(fn, np.array(D)); return D
if __name__ == "__main__":
    import sys; e = int(sys.argv[1]); D = dset(e); print(len(D), D)
