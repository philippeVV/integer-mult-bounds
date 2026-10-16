# TM-native costs of small transposes and masked-lag pass programs: closed as a κ route (×1.00)

Status: scoped closure. A proof for the whole-volume pass model, plus exact negative calibration at small sizes.
Author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`. Multipliers are over κ = 0.0017.

## Verdict

- **Whole-volume masked-lag pass model** ("a pass costs O(V) for any lag"): κ ≡ 0, **proved**. Interchanging two u-bit
  fields needs ≥ log_k(2^{u+1} − 1) passes, coded or not, and routing reaches O(u). The model cannot express the #109 saving
  at any size.
- **Step model** (exact minimum steps of a k-tape oblivious binary TM): XOR never beats routing. This covers 17
  (configuration, size) pairs at step level and 4 at pass level.
- The only gap the frame recursion leaves is leaf fusion, a constant ×0.57–0.60. It multiplies the O(1) term of
  F_k ≤ (s/W)F_{k−1} + O(1), so the exponent is unchanged.
- **Multiplier: ×1.00.**

## Changed hypothesis tested

The subspace-frame contract prices address interchange by s = Σ rank A_e against Wm. Here that price is replaced by
ground-truth machine cost: the exact minimum time of a k-tape binary TM, found by SAT. The question is whether the true
small-size optimum shows a coding (XOR) gain that the frame model misses.

## Proof sketch (lag-set bound, `lagbound.py`)

**Setup.** Let Λ_t be the set of lags ℓ such that some cell c of tape t holds a form whose support contains input c + ℓ.
A pass writing tape o with read lags δ_t gives

  Λ_o′ ⊆ Λ_o ∪ ⋃_{t≠o} (Λ_t + δ_t).

Masks and XOR only combine what was read: cancellation can remove lags but never add them.

**Growth.**
- With one written tape per pass, |Λ_o′| ≤ Σ_t |Λ_t|, so the reachable sizes are k-bonacci.
- With all tapes written per pass, the total grows by at most ×k.

**What the transpose needs.** The 2^{u+1} − 1 lags j(1, −1) on one tape. So passes ≥ log_k(2^{u+1} − 1) = Θ(u/log k), and
the cost per volume grows linearly in u.

**Where #109 escapes.** Its children run on role streams of volume V/W, so lag mass must be counted with volume weights. In
that regime the TM-native question becomes the undirected coding-gap question (AFKL/ZLL). That is the compiled-gadget and
frame-oracle lanes, not a new one.

## Exact tables (scouting run; small cells replayed here)

**Step model.** R×R transpose of b-bit records, workspace L = V + 2. Each entry is route / xor.

| R×R, b | V | k=1, c=1 | k=2, c=0 | k=2, c=1 | k=3, c=0 |
|---|---|---|---|---|---|
| 2×2, 1 | 4 | 4 / 4 | 4 / 4 | 4 / 4 | – |
| 2×2, 2 | 8 | 11 / 11 | 9 / 9 | 9 / 9 | – |
| 2×2, 3 | 12 | 22 / 22 | 15 / 15 | 14 / 14 | 14 / 14 |
| 2×2, 4 | 16 | > 35 (stopped) | 21 / 21 | – | 19 / 19 |
| 2×2, 5 | 20 | – | 27 / 27 | – | 24 / 24 |
| 4×4, 1 | 16 | – | 29 / 29 | 28 / 28 | 25 / 25 |

**Pass model.** 2^u × 2^u interchange.

| u | k | route | xor | lag lower bound | bit-by-bit routing 3u |
|---|---|---|---|---|---|
| 1 | 2 | 3 | 3 | 3 | 3 |
| 1 | 3 | 3 | 3 | 3 | 3 |
| 2 | 2 | 6 | 6 | 5 | 6 |
| 2 | 3 | 5 | 5 | 4 | 6 |

**Leaf fusion.** The 4×4 transpose done as two 1-bit interchanges, each solved exactly:
- k = 2: 25 + 26 = 51 steps against the native 29;
- k = 3: 22 + 20 = 42 against the native 25.

The step model is non-uniform and oblivious with a free program counter. Its optima are lower bounds for real machines, which
is the right direction for a negative calibration.

## Stop test that fired

The test was: stop if XOR never lowers the exact optimum below routing at every solved size, and the pass model is
κ-blind. Both hold.

## Open next question (low value)

The only TM-native variant with any power is volume-weighted: W role streams, with children at cost 1/W. This belongs in the
CP-SAT frame oracle as "row-lag" gates.
- **First experiment:** the 2×2 interchange with W ≤ 6 roles.
- **Stop:** if XOR never beats m = 2. This is expected from expM5: deficits appear only at W ≳ 78.

## Files and reproduction

Run from this directory.

| Script | Run | Runtime | Output |
|---|---|---|---|
| `lagbound.py`: lag-set bound (stdlib) | `python3 -I -B lagbound.py` | 0.1 s | `lagbound.out`, `PASS_LAGBOUND` |
| `replay_small.py`: 9 step-model cells (searched from T = 1, so each minimum is certified by UNSAT below it) and 3 pass-model cells, route = xor | `uv run -q --no-project --python 3.12 --with python-sat python -I -B replay_small.py` | 65 s | `replay_small.out`, `PASS_TM_REPLAY` |
| `tmsat.py`, `passsat.py`: the SAT encodings (python-sat / CaDiCaL), unchanged | — | — | — |

`results-scouting.jsonl` holds the raw records of the full scouting tables (40 rows), including the larger cells that are not
replayed here. Those cells ran 20 s to 2 min each, searched from T_lo = V (4×4: 16 or 26).

Prepared with Claude Code (Claude Opus 5.5) assistance. Not independently reviewed.
