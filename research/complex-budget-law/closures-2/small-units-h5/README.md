# Small cross-field units at h_c = 5 (complement sinks): OPEN, ceiling ×3.03 ('best')

Status: **not closed.** No unit was found, and no certificate of non-existence exists. Author-checked, LOCAL-ONLY until
published. No κ claim. Base: main `0fdd7d82`. Multipliers are over κ = 0.0017.

## Verdict

| Item | Value |
|---|---|
| Best known unit | D6 at x = 1, worth ×0 |
| No unit with x < 1 for v ≤ 24 | **not certified** |
| Ceiling at the LP floor x = 0.4027 | **×3.03** ('best' planning row); ×5.1 on 'cal' |
| Break-even (pair κ = 0.0017) | x ≤ 0.7275 ('best'); x ≤ 0.690 ('cal') |
| ×2 | needs x ≤ 0.567 on either row |
| Multiplier realised | ×1.00 |

So "x < 1" is the wrong target: a unit is worth something only if x ≤ 0.73. The cell is bounded, not closed, and its upside
is at most ×3.

## Setting and hypothesis

This is the last open cell of the small-unit lever, in the complex orientation at h_c = 5:
- odd labels q_T and odd complement sinks w_T in F2^5, with w_T·q_T = 1 (136 types);
- a decoder N = Σ_j r_j a_jᵀ over Q, asymmetric and not Gram, with N∘B = diag;
- totals of rank β_j = rank_F2 B[R_j, A_j], which may be below full.

The question: is there a unit with x = Σβ/v < 1 for v ≤ 24?

## What was tried (none conclusive)

**Type-level integer model (CP-SAT, then HiGHS with symmetry-closed rows).**
- The model uses the H (Hadamard) and A (acyclic) rows. It never became infeasible.
- Cheap integer points keep appearing, e.g. a single β = 5 total over 12–16 types.
- Closing the rows under the 1,440-element symmetry group runs out of memory at the 4G cap.
- The H and A rows are scale-free, with fractional optimum 0.4027. A certificate needs realisability rows: per-port
  cancellation, minrank.

**Port-level CP-SAT search for actual units.**
- Every run is UNKNOWN: v = 6 (120 s); v = 8, 12, 16, 20 with coefficients in {−1, 0, 1} (300 s each); and the D6-pairing
  puncture search (300 s).
- Bounded-coefficient results are not certificates in any case, because Q decoders are unbounded.
- **Sanity check:** with x ≤ 1 allowed, the solver finds an exact unit in under 1 s (`hit_v6_k3_c1.json`: v = 6, three β = 2
  totals, x = 1, h_c = 4, h_b = 3). `verify_hit.py` re-verifies it exactly.

## Related cells

- Standard-sink cells at h ≤ 6 and complement-sink cells at h ≤ 4 are closed for ×10 (expN2, expN4).
- Complement sinks at h_c = 6 are bounded only by the ratio bound x ≥ 4/15 (expN5). That cell was not tested and is also
  open.

## Next step (only if a formal close is wanted)

1. **Certify against break-even (x ≤ 0.7275), not x < 1.**
   - Method: port-level lazy rows (pair cancellation: 0 or ≥ 2 covers), plus H and A restricted to the touched ports, with
     the types fixed per orbit representative.
   - **Stop:** UNSAT for v ≤ 24 on two solvers. Then the cell is closed for any gain.
2. **Otherwise, test whether the ILP's cheapest integer point is realisable** by exact linear algebra.
   - **Stop:** drop the cell if it is not realisable with integer coefficients up to 3.

## Files and reproduction

Run each from this directory: `python3 -I -B <script>`. Everything is stdlib.

| Script | Runtime | Output |
|---|---|---|
| `score.py`: pair scoring at h_c = 5 (bit twin h_b = 4, same x), 'cal' and 'best' rows, and the break-even values. Uses `pair_s3.py` and `transfer.py` (expNE) and `kappa_pair.py` (expNB), copied with only their import paths changed. | 3.8 s | `score.out`, `PASS_EXPNF_SCORE` |
| `verify_hit.py`: exact verifier of a unit and its bit twin (uses `unit.py` from expNA), run on the sanity hit | 0.2 s | `verify_hit.out`, `PASS_EXPNF_VERIFY` |

Not included: the integer-model scripts (`model5.py`, `group5.py`, `ilp.py`, `ilp2.py`) and their caches (≈ 53 MB), and the
CP-SAT unit searches (`cpunit.py`, `cpdec.py`). Their results are all UNKNOWN or out-of-memory, and they are quoted above.

Prepared with Claude Code (Claude Opus 5.5) assistance. Not independently reviewed.
