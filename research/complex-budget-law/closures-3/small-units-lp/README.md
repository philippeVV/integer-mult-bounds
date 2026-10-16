# Small units, LP rows (update to closures-2/small-units-h5): trace rows are invalid; the lane is PAUSED

Status: update to `../../closures-2/small-units-h5/` (still **OPEN**, now **paused**). Author-checked, LOCAL-ONLY until
published. No κ claim. Base: main `0fdd7d82`. Multipliers are over κ = 0.0017.

## Verdict

| Item | Value |
|---|---|
| Is "tr X_j ≤ ρ′_K" (or "v ≤ Σ_j min(k_j, ρ′_j)") a valid general LP row? | **No.** Counterexample (`spread.py`, exact): with one unit of slack, every total is non-idempotent and Σ rank X_j = v + 1. Per-total idempotency is worth only **O(1)** ports per unit unless a pattern lemma forbids defect spreading (open) |
| Where the per-total row is valid | only for **rank-additive** units (Σ rank X_j = v, slack 0) |
| Exact idempotent-rank caps ρ′_K | h = 6: K* (orbit 26) and orbit 32: 8; Fano β = 3 orbits: 4; **orbit 77 (β = 5, 16 × 32): ρ′ = 16 = rank_Q** (exact witness, `idem77.py`). h = 5: full β = 5 block 15; Fano 4 (`rpub.py 5 7 4`) |
| Tight-regime LP (rank-additive units only) | h = 6: x ≥ 5/16 = 0.3125, and **no idempotent-cap row can certify more**, because orbit 77 is an exact rank-16 idempotent. h = 5: x ≥ 1165/2837 = 0.41065 (up from 1085/2669) |
| Pair rows rank_Q N ≤ h_b | **cannot move the LP** (absolute caps on a homogeneous LP); they bite only with integrality (v ≤ 17h_b at h = 6, ≤ 16h_b at h = 5). K* and the h = 5 full block stay coverable at h_b = 5 |
| Stop tests (x ≥ 0.45 at h = 6; ≥ 0.49–0.50 at h = 5) | **both FAIL** |
| Unconditional floors | unchanged: x ≥ 4/15 (h = 6) and x ≥ 1085/2669 (h = 5) |
| Multiplier | realised ×0 (no unit with x < 1 is known); planning ceilings ×16.7 (h = 6, 5/16) and ×14.7 (h = 5, 1165/2837) for rank-additive units, against the closures-2 pair-scored ×3.03 ('best') at h_c = 5 |

**Decision: pause the small-unit LP lane.** The cheap certificate routes (trace rows, idempotent caps, pair rows) are exhausted
or invalid; what remains needs either an orthogonality test on orbit-77 idempotents or an integer program in the pair regime
(both below). Reopen if either of those reports, or if a unit with x ≤ 0.7275 (the closures-2 break-even) is exhibited.

## Changed hypothesis

expN9's next step 1: add idempotency (trace) rows to the small-unit LP, plus the pair constraint rank_Q N ≤ h_b. The hope was
that "each total is an idempotent of bounded rank" would lift the LP floor past the stop lines.

## The idempotency row: what is and is not valid

Write each total X_j = U_jV_jᵀ (rank k_j), U = [U_1 … U_m], V likewise, so UVᵀ = I_v; K = Σk_j.
- E := VᵀU is a K × K idempotent of rank v, and F := I − E has rank s = K − v (the H-row slack).
- X_j² − X_j = U_j(E_jj − I)V_jᵀ, so **rank(X_j² − X_j) ≤ s** for every j, and Σ_j tr X_j = v.
- **Valid:** v ≤ Σ_j k_j − max_j δ(k_j), where δ_K(k) ≥ 1 iff k > ρ′_K. A gain of O(1), not Θ(v).
- **Invalid:** "tr X_j ≤ ρ′_K" (traces of non-idempotent totals are unbounded) and "v ≤ Σ_j min(k_j, ρ′_j)". In `spread.py`,
  K = v + 1, so F = fgᵀ has rank 1, every block has rank(X_j² − X_j) = 1, and the defect spreads over all totals at total cost 1.
- The LP rows are homogeneous (v → ∞), so no per-total idempotency row is valid there unless the D_r B̃ D_a patterns forbid
  spreading.

## Idempotent-rank caps (exact)

X is idempotent iff B̃′ = B̃′LB̃′ with L on type positions, i.e. I_k ∈ span{R[:, q]C[w, :]ᵀ} for a rank factorisation B̃′ = CR.
Feasibility is monotone under adding sinks or labels that do not raise the rank, so ρ′ ≤ c needs only the rank-maximal supports
of rank > c (`rpub.py`, mod p first, then exact over Q with FLINT).

| h | orbit | β | block | rank_Q | ρ′ | how |
|---|---|---|---|---|---|---|
| 6 | 26 (K*) | 4 | 16 × 16 | 15 | **8** | `rpub.py 6 26 8`: 17,016 rank-maximal supports, none feasible (about 27 min, not replayed) |
| 6 | 32 | 4 | 16 × 16 | 15 | **8** | `rpub.py 6 32 8` (not replayed) |
| 6 | 77 | 5 | 16 × 32 | 16 | **16** | full support feasible; `idem77.py`: X on 256 ports with λ ∈ {1/8, −7/8}, X² = X, rank 16 |
| 6 | 16, 18, 20 (Fano) | 3 | 8 × 16 | 7 | **4** | `rpub.py 6 k 4` (not replayed) |
| 5 | 43 (full) | 5 | 16 × 16 | 16 | **15** | the rank-16 support is the whole block, infeasible (expN8) |
| 5 | 7 (Fano) | 3 | 8 × 8 | 7 | **4** | `rpub.py 5 7 4` (replayed) |

The Fano gain (7/3) and the K* gain (15/4) are not idempotent-compatible (they fall to 4/3 and 2), but the β = 5 block at h = 6 is
idempotent at full rank 16, and at h = 5 the full block keeps 15/5 = 3. `capsdef.py` holds the cap table used by the LP.

## Tight-regime LPs (I row added; rational duals; scouting results, not replayed)

- **h = 6:** certificate x ≥ 5/16 from the single I row. The primal moves from K* to orbit 77 (m uniform on the 480 cross types,
  y = 1/16 on the 64 β = 5 totals). With caps at the lower bounds the value is the same, so 5/16 is the most an idempotent-cap
  row can certify at h = 6.
- **h = 5:** LP 0.410645; exact H separation and CP-SAT A separation found no violated row; rational dual x ≥ 1165/2837 (I row
  at weight 710/2837 and 6 A rows).

| h | bound | κ ≤ | × over 0.0017 | stop line |
|---|---|---|---|---|
| 6 | 5/16 (rank-additive) | 2.83e-2 | ×16.7 | 0.45 |
| 5 | 1165/2837 (rank-additive) | 2.50e-2 | ×14.7 | 0.494 / 0.50 |

## Next steps if the lane is reopened

1. **Orthogonality of orbit-77 idempotents** (h = 6): solve exactly over Q for λ's with X_i² = X_i and X_iX_j = 0 across the 64
   orbit-77 blocks (bilinear; fix one block to the witness and solve the rest linearly). **Stop:** if infeasible for every pair,
   rerun the I-row LP with orbit-77 pairs excluded; if that value is > 0.45, h = 6 is closed for rank-additive units.
2. **Integer program in the pair regime** (h = 5, h_b = 5: v ≤ 80, multiplicities ≤ 5, α(S) ≤ 5, integer totals, slack
   s ≥ max_j δ(k_j)). **Stop:** integer optimum > 0.50 closes h = 5 at ≤ ×2.4 realistic.
3. Do not pursue trace rows in the homogeneous LP further.

## Files and reproduction

Run each from this directory. Runtimes under a 4G cap, one worker.

| Script | Checks | Runtime | Output |
|---|---|---|---|
| `uv run -q --no-project --python 3.12 --with python-flint python -I -B spread.py` | the slack-spreading counterexample (exact): Σ rank X_j = v + 1, every block non-idempotent | 0.1 s | `spread.out`, `PASS_EXPNG_SPREAD` |
| `uv run -q --no-project --python 3.12 --with numpy --with python-flint python -I -B idem77.py` | exact rank-16 idempotent on orbit 77 (h = 6) | 96 s (mostly rebuilding the h = 6 orbit data) | `idem77.out`, `PASS_EXPNG_IDEM77` |
| `uv run -q --no-project --python 3.12 --with numpy --with python-flint python -I -B rpub.py 5 7 4` | exact cap ρ′ ≤ 4 on the h = 5 Fano orbit | 1.4 s | `rpub_5_7_4.out`, `PASS_EXPNG_RPUB` |
| `core6.py`, `rhoprime.py`, `capsdef.py` | closed columns and orbits (now rebuilt instead of cached to a pickle), the exact idempotent-support test, the cap table | — | — |

Not included: the LP drivers and certificates (`symlp.py`, `symlp5.py`, `icap5.py`, `cert5.py`), the lower-bound searches
(`rpgrow.py`, `lball.py`), the pair coverability CP-SAT (`pair.py`, `pair_k.py`) and their caches (about 22 MB). Their results are
quoted above. The h = 6 `rpub.py` scans run the same way (`rpub.py 6 26 8`, about 27 min).

Source: expNG. Prepared with Claude Code (Claude Opus 5.5) assistance. Not independently reviewed.
