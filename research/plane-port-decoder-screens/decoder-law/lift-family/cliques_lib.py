from e8 import *
NO = [sum(1 << b for b in range(NL) if b != a and IP[a][b] != 0) for a in range(NL)]
def random_clique(k, rnd):
    while True:
        S = [rnd.randrange(NL)]; cand = NO[S[0]]
        while len(S) < k and cand:
            c = rnd.choice(list(bits(cand))); S.append(c); cand &= NO[c]
        if len(S) == k: return S
def common(S):
    """C(S) = lines b with b in S or b non-orth to every s in S (the allowed B)."""
    return [b for b in range(NL) if all(b == s or IP[b][s] != 0 for s in S)]
