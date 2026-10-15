# QQ16 final contract: ordered decoders for banks 0–3 and an actual-cover bank-1 loss bound

Status: **result package, local draft; nothing pushed or posted.**
- Card `joint-query-controller`, priority revision 4, base main `0fdd7d82`.
- Follow-up to the final-record census (fork branch `findings/joint-query-first-pair`,
  `research/joint-query-first-pair-catalogue/final-record/`, #2074 comment 18950410).
- That census closed with two open items: "no free-algebra decoder" and "no actual-cover witness".
  This package addresses both.

**Contract (PR2346).** FINAL is:
- the final selected X/Y fields of all four pairs;
- the final Y ordinary fields of all pairs (completed-X descendants);
- the forward histories f0, f2, f4, f6, f8.

INIT64 is B0, B2 applied to X_0U, and B_(2k+1), B_(2k+2) applied to X_kU for k = 1, 2, 3.
d_A,k is Y_kA just after bank k's first backward call (events 1, 5, 9, 13).

The maps are CONSERVED-QQ16 (1)/(3). The bridges and terminal follow PR2346 §1 and PR2174 §2.
The schedule is the pinned 19-event recipe, with local pairs 0,0,0,1,1,2,2,3,3.

## Results

### 1. Explicit ordered decoder for bank 0: FINAL + INIT64 + d_A,0 (task 1, done)

This is PR2174's decoder (5) with its cut inputs replaced by final-contract quantities. Three
facts justify the replacement:
- **Conservation.** F0, B1 and F2 preserve X_p+Y_p on A∪T, and bridges and the terminal touch
  U only. So the cut sums are p_A = X_0A^fin + Y_0A^fin and p_T = X_0T^fin + Y_0T^fin.
- **Tail identity (PR2174 §2).** It gives C_0 = Y_0U^fin = w_0 + K_2 f_2.
- **The stated INIT64 reads.** q_0 = B_0 X_0U and q_2 = B_2 X_0U.

```text
w0 = C0 - K2 f2
t  = pA - dA0 + P_2A [ B2 w0 + q2 - B2 K0 f0 + ell_2T (pT - P_0T f0) ]
xT = P_0T ell_0A t
xA = P_0A [ f0 - q0 - ell_0T xT ]
yA = pA - xA,   yT = pT - xT
```

- Every product is ordered.
- There is no inverse of an ordinary query, route, centre response or Schur complement.
- Only the dyadic sections are used.
- f4, f6, f8, the final fields of the other banks and the X_0A^fin/X_0T^fin split are not needed.

### 2. Banks 1–3 (task 3, done)

The tail identity holds for every bank, for arbitrary inputs in all banks:

    C_1 = w_1 + K4 f4,   C_2 = w_2 + K6 f6,   C_3 = w_3 + K8 f8.

The Y_kU that enters bank k's first backward call is V_k + X_kU^init, where
V_k = w_k − C_(k−1) + C_(k−2) − … . So B_s and B_t of it are known from INIT64.

For bank k, with s = 2k+1 and t = s+1:

```text
beta = ell_sT P_tT ell_tA X_kA^fin                  (= b + 12a of B_s; X_kA^fin = P_tA ell_tT P_sT beta)
Bs_w = B_s V_k + INIT(k,s),   Bt_w = B_t V_k + INIT(k,t)
xT = P_sT [ (beta - Bs_w - ell_sA dA_k - ell_sT pT) / 11 ]
yT = pT - xT
yA = dA_k - P_tA [ Bt_w + ell_tT yT ]
xA = pA - yA
```

- **Banks 1 and 2** need d_A,1 and d_A,2. Their f4 and f6 are useless for this purpose: the
  bridge before F_t adds X_(k+1)U, so f_t contains B_t X_(k+1)U, which no INIT64 read sees.
- **Bank 3 needs no record.** No bridge separates B7 from F8, so
  `dA_3 = pA + P_8A [ INIT(3,8) - Bt_w + ell_8A X_3A^fin - f8 ]`.
- **The 1/11 is inherent to this contract.** x_T reaches the observations only as
  12a − a = 11a inside beta. It is an odd denominator (precision note). PR2174's bank-0 decoder
  has none.

So **FINAL + INIT64 + d_A on banks 0, 1, 2 (88 arrays) determines all 128 selected inputs.** The
decoder is explicit, which turns the census's sufficiency claim into a derivation.

**Equivalent all-initial contract.** Replace the three d_A records by three initial cross-bank
reads: B2 X_1U, B4 X_2U and B6 X_3U. f_t then yields d_A,k, for example
`dA_0 = pA + X_0A^fin − P_2A[f2 + B2 w0 − B2 X_1U]`. This is again 88 arrays, and no mid-run
selected tap is needed. Under storage-only pricing it is not cheaper, because the reads are live
from event 0.

### 3. Actual-cover witness for the bank-1 loss under FINAL + INIT64 + d_A,0 (task 2, done)

**Exact reduction.** Every z ∈ ker B5 ∩ ker B6, applied to X_2U, gives a full input as follows.
Put X_2U = 11z and u = 11 P_4A B4 z. Set bank 1's selected inputs to:

```text
xT = -P_3T ell_3A P_4A B4 z,  yT = -xT,  yA = u - P_4A ell_4T yT,  xA = -yA
```

Set everything else to zero. Then every observation in FINAL, INIT64 and d_A,0 vanishes.

Trajectory:
1. β = 0 after B3.
2. f4 = B4(11z) − 11 ell_4A P_4A B4 z = 0.
3. Bank 1's ordinary residue K_3U a sits in Y_1U. It leaves only through final X fields, which
   are unobserved.
4. Y_2U = 11z enters B5, which sees only B5 z = B6 z = 0.

Bank 1's selected input is nonzero whenever B4 z ≠ 0, and d_A,1 = u then sees it. Bank 2 is the
same with (6; 7, 8).

**Actual routes.** `routes_bound.py` reads the actual chart routes D_v(p) = M(g, δ), taken from
the first 112 X-alignment cover maps of each F_v program in the pinned certificate:
- all are involutions;
- their owner parts match the Y alignments of every B_s_t program (1792 maps);
- for each ordered pair among charts {4, 5, 6} (and {6, 7, 8}, {2, 3, 4}), the 96 relative
  owner routes D_a(p) D_b(p)^−1 are pairwise distinct.

Then G = [B4;B5;B6][B4;B5;B6]* has tr G = 1152 and tr G² = 64512, exactly. Cover maps act freely,
so every nonidentity route word has zero trace, and Cauchy–Schwarz gives
rank[B4;B5;B6] ≥ 1152²/64512 = 144/7 > 16 ≥ rank[B5;B6]. Hence

    dim (bank-1 selected inputs invisible to FINAL + INIT64 + d_A,0) ≥ dim B4(ker B5 ∩ ker B6) ≥ 32|G|/7,

in PR2346's normalization (per owner: ≥ 32/7 of 8, against exactly 8 in generic matrix models).
Bank 2 has the same bound. d_A,1 is 8 arrays and repairs the loss, so the loss lies between 32/7
and 8 arrays per owner.

- This uses the same two-moment method as PR2346 (Ross Richardson #2104), applied to a
  different triple.
- PR2346's bank-0 bound, ⌈8|G|/7⌉ (an intersection of the (1;0,2) and (2;3,4) images), is
  reproduced as a control. Banks 1–2 need no intersection, so their bound is 4× larger.
- This is an existence proof. It gives no compact witness vector, as in PR2346.

### 4. Correction to the final-record README, item 4

The "discrepancy" there is an artifact.
- PR2174's own contract recovers bank 0 in both models (INIT B0/B2 on X_0U, cut sums p_A and
  p_T, d_A, f0, f2, C_0). Dropping C_0 loses 8.
- The census rows labelled pA0/pT0 in `modelq2.py` are the cut X_A/X_T, not the sums.

Item 4 should be withdrawn.

## Evidence (check.py → `PASS_QQ16_FINAL_CONTRACT_DECODER`, about 3 min)

| part | what | status |
|---|---|---|
| A | all 7 decoders (banks 0–3, plus the all-initial variant) recover every selected input exactly. Random inputs on all four banks, two primes (2^31−1 and a 62-bit prime). Routes: random noncommuting m×m matrices (m = 3, 8), D_0 = I and D_0 random, plus permutation, diagonal and identity families. 24 runs. | randomized exact substitution |
| B | 30 decoder mutations rejected: tail sign or chart, dropped K0 f0, wrong section chart, dropped alternating C, dropped INIT, 11→12, wrong d_A. An order-swapped β passes only with commuting routes, so the test sees operator order. | control |
| C | simulator mutations detected: bridge sign, pair = stage//2 mod 4, backward 11→12, 12→13 | control |
| D | rank census. Generic matrix model, m = 2: lost per owner is [8,8,8,0] → +d_A,0: [0,8,8,0] → +d_A,0..2: 0; the +XB contract: 0. Permutation model, m = 3: reproduces modelq2's 16/3 exactly, from an independent simulator. Both primes agree. | finite model |
| E | witnesses assembled from a computed kernel vector and fully simulated (matrix and permutation models, two primes). Banks 1 and 2 are blind to FINAL+INIT64+d_A,0. The PR2346-type bank-0 witness is blind to FINAL+INIT64. Each is seen by its own d_A,k. | finite model |
| F | actual-route bound: exact rational traces from the pinned certificate's 8×112 routes, with a D5 := D4 control that collapses to 13.5 ≤ 16 | **actual cover** (premises below) |

**What is and is not established.**
- **Decoders.** The decoder identities are short hand derivations, written above.
- **Decoder checks.** They are checked by randomized exact substitution with noncommuting routes.
  They are not proved symbolically or in a proof assistant.
- **Actual-cover bound.** It is finite linear algebra, but it inherits two premises from PR2346:
  - the logical route model, in which D_v(p) is the owner/translation cover map (the native
    phase frames are abstracted, as in PR2174 and PR2346);
  - free action of the cover maps.
- **Physical cost.** Nothing about physical access, storage or cleanup of the d_A records (or of
  the XB reads) is priced here. The price gate's convention question (rank-8 copy versus one
  unit per retained array per call) is unchanged.
- **Review.** A sub-agent did the work. It has not been independently reviewed.

## Run

```sh
export REPO=/path/to/main/checkout     # needs research/joint-query and the selected-replay archive
uv run --python 3.12 --with python-flint python check.py          # ~3 min (D census dominates)
uv run --python 3.12 --with python-flint python check.py --quick  # ~4 s, skips D
uv run --python 3.12 --with python-flint python routes_bound.py   # ~3 s, part F alone
```

Files:
- `qq16mat.py`: the operator model and simulator (no predecessor executable imported).
- `decode.py`: the decoders.
- `check.py`: the checker.
- `routes_bound.py`: the actual-route bound.
- `CHECK.json`, `check.log`: the last run.

## Credit

- Ali Romman (#1969): the QQ16 maps.
- commdata2338 (PR2174): the decoder (5) and the tail identity, extended here to all banks.
- alephnan (PR2346): the contract, and the actual-cover theorem whose method this follows.
- Ross Richardson (#2104 two-moment method; PR2480).
- thakrarsagar (PR2470).
- ywz65 (PR2468/2492).

Prepared by Philippe Van Velzen with substantial Claude Code (Claude Opus 5.5) assistance. Apache-2.0.
