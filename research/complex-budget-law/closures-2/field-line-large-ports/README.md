# Field-line (polar-space) large ports: stop test fired (best exact γ = 1, ×1.00)

Status: scoped closure for field-line label systems. The decoders are exact; the residual minrank law is **numerical
evidence only**. Author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`. Multipliers are over
κ = 0.0017.

## Verdict

| Quantity | Value |
|---|---|
| Best exact decoder, all 8 families | **γ = 1** (`exact.py`), so ν ≤ 0 and there is no gain |
| Multiplier | **×1.00** |
| Unexcluded optimistic ceiling at the proven floors | ×24–36 (`kappa.py`) |

The gap between floor and achievable is explained by a numerical law (below), which is **not proved**.

## Changed hypothesis tested

Ports are u-dimensional, with label fraction k = u/(h−1) ≈ 1/4–3/8. Each label is an F_q-line (q = 2^u) of a classical polar
geometry, read as a u-dimensional F2-subspace under Tr∘B. The families:
- PG(n−1, q) with the dot or Hermitian form;
- the E8-style "+1" chart AG(n−1, q).

The pairing ranks are in {0, u}, so L1′ applies: C ≥ u·minrank ≥ u·α.

## Exact decoders (`exact.py`)

- **D0 (identity):** γ = 1.
- **D1 (frame decoder):** a CP-SAT maximum packing of pairwise-free ω-sets (self-polar frames), N = Σ 1_C 1_Cᵀ. It is exactly
  valid and realises minrank about V/ω. Every frame spans F_q^n, so each frame term costs h, and **C = uV exactly: γ = 1** in
  all 8 families. This is the u-dimensional form of PR #2244 Theorem C (unbiased frames give x = 1).
- The best invariant PSD decoder (Hoffman) has rank of about V/2 or more, so γ_dense = 1.6–2.2.

| Family | h | u | V | ω | α floor (γ ≥) | Best exact γ |
|---|---|---|---|---|---|---|
| PG(2,4) dot | 6 | 2 | 16 | 3 | .375 | 1 |
| PG(3,4) Hermitian | 8 | 2 | 40 | 4 | .175 | 1 |
| PG(3,4) dot | 8 | 2 | 64 | 4 | .156 | 1 |
| AG(3,4) dot +1 | 8 | 2 | 48 | 4 | .208 | 1 |
| AG(2,8) dot +1 | 9 | 3 | 56 | 3 | .286 | 1 |
| PG(2,16) Hermitian | 12 | 4 | 208 | 3 | ≥ .212 | 1 |
| PG(2,16) dot | 12 | 4 | 256 | 3 | ≥ .203 | 1 |
| PG(3,8) dot | 12 | 3 | 512 | 4 | ≥ .078 | 1 |

## Residual law: NUMERICAL, not proved

**The law: "masked minrank = V/ω" for polarity graphs.** In every family where Levenberg–Marquardt was run below the frame
value, the best residual is exactly V − ωK, and the fit converges only at K ≈ V/ω. Examples:
- PG(3,4) dot: residuals 32, 24, 16, 8 at K = 8, 10, 12, 14;
- PG(3,4) Hermitian: residuals 12, 8, 4 at K = 7, 8, 9.

**Why it is plausible.** The free graph of field lines has density about 2^{−u}, while E8's has about ½. E8's super-frame
compression (V/K = 14 against ω = 8) has no analogue here.

**Status.** These are floating-point fits, uncertified. If the law holds, every decoder has γ ≥ u/h ≈ 0.25–0.33. That alone
sits at or above the ×10-calibrated thresholds (γ ≤ 0.21 at k = 1/4; 0.31–0.34 at k = 3/8).

## Small exact lemma (Bessel bound)

**Statement.** A PSD decoder whose support graph is triangle-free has rank ≥ V/2.

**Proof.** Each row has Σ G_ab² ≤ 2, so V²/rank ≤ ‖G‖_F² ≤ 2V.

**Consequence.** For field-line labels with h < 3u, every PSD decoder has γ ≥ ½.

## Stop test that fired

The test was: stop if the best exact decoder has γ > 0.705 in every family. It fired (γ = 1 everywhere).

The second condition, layered floors ≥ 0.51, did not fire: the floors are 0.16–0.29.

## Open next questions

1. **Prove the minrank law** (e.g. a Haemers/ϑ-type bound using the 2^{−u} free density). That would close k ≥ 1/4 for all
   field alphabets. **Stop:** if nothing beyond α/V comes out in one session.
2. **Non-field u-dimensional labels** (pairing ranks in {1, …, u−1}).
   - Prerequisite: free-pair density ≳ ¼ with ω ≤ h/u.
   - **Stop:** if no structured system at h ≤ 17, u ≥ 3 has density ≥ 0.1.

## Files and reproduction

| Script | Run | Runtime | Output |
|---|---|---|---|
| `exact.py`: exact D0/D1 decoders with F2 pairing ranks (numpy, ortools CP-SAT with a 120 s limit) | `uv run -q --no-project --python 3.12 --with numpy --with ortools python -I -B exact.py <u> <n> dot\|herm [plus1]` for the 8 families: `2 3 dot`, `2 4 herm`, `2 4 dot`, `2 4 dot plus1`, `3 3 dot plus1`, `4 3 herm`, `4 3 dot`, `3 4 dot` | about 1 s each for the first five; 121–123 s each for the last three (CP-SAT time limit) | `exact.out` (all 8 concatenated; 8 × `PASS_EXPNM_EXACT`) |
| `kappa.py`: ceilings via `../../framework-ceiling/framework_ceiling.py` (stdlib) | `python3 -I -B kappa.py` | 0.1 s | `kappa.out`, `PASS_EXPNM_KAPPA` |
| `fam.py`: GF(2^u), PG/AG points, forms, F2 labels | — | — | — |

The LM screen (`screen.py`) and its long runs are not included. Their results are quoted above as numerical evidence.

Prepared with Claude Code (Claude Opus 5.5) assistance. Not independently reviewed.
