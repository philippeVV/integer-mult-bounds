#!/usr/bin/env python3
"""Step 1 of the bit plane-port experiment: exact checks behind the class-P theorem and the general rank certificate.
Standard library only.  usage: python3 -I -B verify_step1.py [--full]
  default: 11,950 (|S|=5, rank M=4) sets for F1/F2 + clique certificates; census/cost-lemma on samples; s=4 small-S check on a sample.
  --full : F1/F2 on all 226,903 maximal sets, census on all 11,950, cost lemma on 2 sets, small-S on all 3- and 4-sets."""
import sys, json, itertools, time
from collections import Counter
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import core
FULL = '--full' in sys.argv
t0 = time.time()
L = list(range(111)); P = list(range(111, 120))
# --- A. conventions: incumbent P111 decoder costs 63 under rational pairing rank (Theorem C, bit side) ---
E = core.line_E(); Lset = set(L); tot = 0; f2tot = 0
for k in range(9):
    y = [int(i == k) for i in range(9)]; z = core.QF[k]
    R = sorted(core.codeword(y) & Lset); A = sorted(core.codeword(z) & Lset)
    tot += core.pair_rank(R, A, E)
    f2tot += core.rank2([sum(core.B0[r][a] << j for j, a in enumerate(A)) for r in R])
assert tot == 63 and f2tot != 63
print(f'A  incumbent cost: rational pairing rank 63 (F2 label-pairing rank would give {f2tot}); B0 rank_F2', core.rank2([sum(core.B0[i][j] << j for j in range(120)) for i in range(120)]))
# --- B. enumerate pairwise-odd common-r maximal sets (as probe_plane.py) ---
import sets as _sets_mod  # noqa: writes sets_all.json; prints class counts
allsets = json.load(open('sets_all.json'))
fam = [s for s in allsets if len(s[1]) == 5 and s[2] == 4]
assert len(fam) == 11950 and len(allsets) == 226903
# --- C. F1/F2 premises ---
import premises as PR
target = allsets if FULL else fam
seen = {}; bad = 0
for r, S, rk in target:
    key = tuple(S)
    if key not in seen: seen[key] = PR.F1(S)
    bad += (not seen[key]) + (not PR.F2(r, S))
assert bad == 0
print(f'C  F1 (X forced = QF on L minus S) and F2 (Y_a = Z_a = 0) hold on {len(target)} sets')
# --- D. Lemma-1 census on L minus S: reduced reader spans are exactly lambda_x, x in P u S (no D7, dims >= 7) ---
import reduced as RD
samp = fam if FULL else fam[::600]
for r, S, rk in samp:
    mind, nred, ok = RD.census(S); assert mind == 7 and nred == 14 and ok
print(f'D  reduced-side census exact on {len(samp)} sets: min dim 7, reduced sides = 14 root functionals of P u S')
# --- E. clique (F2 rank) certificates for the plane target, any supports ---
import cliqueall  # noqa: asserts every pair of Z non-orthogonal; prints distribution
print('E  (above) plane-masked rank >= 12 certificate for all 11,950 sets; line control omega 8')
# --- F. small-S residual check (|S| = 4: ten cost-6 label terms impossible) ---
import small_s as SS
S4 = sorted({c for r, S, rk in allsets for c in itertools.combinations(S, 4)})
S4 = S4 if FULL else S4[::200]
cnt = Counter(SS.check_set(S, 4)[0] for S in S4); assert set(cnt) <= {'nosol', 'none'}
print(f'F  |S|=4: no ten-term cost-6 label decoder on P u S for {len(S4)} sets {dict(cnt)}')
if FULL:
    S3 = sorted({c for r, S, rk in allsets for c in itertools.combinations(S, 3)})
    cnt = Counter(SS.check_set(S, 3)[0] for S in S3); assert set(cnt) <= {'nosol', 'none'}
    print(f'F  |S|=3: no ten-term label decoder of cost <= 61 on P u S for {len(S3)} sets {dict(cnt)}')
import small12 as S12
S2 = sorted({c for r, S, rk in allsets for c in itertools.combinations(S, 2)})
S1 = sorted({c for r, S, rk in allsets for c in itertools.combinations(S, 1)})
if not FULL: S2 = S2[::300]; S1 = S1[::20]
for s_, lst in ((2, S2), (1, S1)):
    cnt = Counter(S12.check(S, s_)[0] for S in lst); assert set(cnt) == {"none"}
    print(f"G  |S|={s_}: no ten-term label decoder with the required savings on P u S for {len(lst)} sets")

# --- H. class U (one side of a term inside S may be arbitrary): F3 for r live; r punctured needs the n1 = 9 / 10 checks ---
lab = [sum(core.U[t][i] << i for i in range(9)) for t in range(120)]
f3 = Counter()
for r, S, rk in (allsets if FULL else fam + allsets[::50]):
    vs = [lab[b] for b in L if b not in set(S) and core.G[b][r] != 0]
    f3[(r >= 111, core.rank2(vs))] += 1
assert all(k == (False, 9) or k == (True, 8) for k in f3)
print(f'H  Odd_r label rank: {dict(f3)}  (9 for every live r; 8 = kernel lambda_r for punctured r)')
import sav9fast, ten19, control19  # noqa: control19 builds an explicit savings-19 ten-term decoder and checks ten19 finds it
assert sav9fast.max_savings(P) == 9
samp = fam if FULL else fam[::400]
worst = max(sav9fast.max_savings(P + S) for r, S, rk in samp)
assert worst <= 12
no19 = all(not ten19.exists19(P + S)[0] for r, S, rk in samp)
assert no19
print(f'H  r punctured, |S|=5: max 9-term savings {worst} <= 12 and no 10-term savings-19 configuration on {len(samp)} sets')
print(f'PASS_BIT_PLANE_STEP1_CLASSP ({time.time() - t0:.0f} s)')
