# Lemma S: support SAT for co-star decoders with k = 7–9 normals (sampled cases)

**Scope.** Exact UNSAT results for the *sampled* cases below. The complete k = 7 family has since
been enumerated and closed in `../k7/`, which uses this lemma, with capacity 3 allowed, on all 2,342 k = 7
engine survivors. The k = 8 and 9 families have **not** been enumerated, so k = 8 and k = 9 remain
open as theorems (except the sets containing F0, `../f0/`). No κ, decoder or cost claim.

## Lemma S

Model and notation as in `../README.md`: normals a_g, capacities c_g ∈ {1,2} with Σ c_g = 9, and
group g contributes `D_g T_gᵀ` with supports in the co-star `S_g`. Introduce Booleans

- `A[g,r] ⇔ D_g[r,:] ≠ 0` and `B[g,r] ⇔ T_g[r,:] ≠ 0` for r ∈ S_g. These are invariant under
  GL(c_g) inside a group.

Every decoder satisfies:

- **(D)** For each live r, some g has `A[g,r] ∧ B[g,r]`, since `M[r,r] = 4`.
- **(X)** For r ⊥ s and a capacity-1 group x with `A[x,r] ∧ B[x,s]`, some other group h has
  `A[h,r] ∧ B[h,s]`. Otherwise `M[r,s]` would be a single nonzero product.
- **(R)** For an orthogonal clique Q, `|Q| = rank 4I_Q ≤ Σ e_g` over groups touching Q (some A
  and some B on Q), where `e_g = min(c_g, |S_g ∩ Q|)`. With `--frames`, only the 1,260 live
  8-frames are used. By the Bessel bound `|S_g ∩ Q| ≥ 4`, so e_g = c_g and the slack is 1. Every
  capacity-2 group must then touch Q, and at most one capacity-1 group may fail to.

Auxiliary variables are defined in one direction only, so every real decoder induces a satisfying
assignment. UNSAT therefore proves that no decoder exists for that (normals, capacities) case,
over any field of characteristic ≠ 2. No Nullstellensatz certificate and no triangle lemma are
needed.

## Results

`check_S.py` is standalone. It rebuilds E8 by root enumeration, pins #1320's line-order SHA-256,
enumerates all 262,928 orthogonal cliques with its own bitmask search, and solves every case
with CaDiCaL 1.5.3 and Glucose 4.

| cases | source | result |
|---|---|---|
| `cases_k7.json` (53) | every k=7 engine survivor found: the 7 from the previous addendum plus 46 from two new uniform samples of 4,000 draws each (≈70k engine pairs) | all UNSAT on both solvers, `--frames` |
| `cases_k9a.json` + `cases_k9b.json` (416) | every engine survivor from 1,030 sampled k=9 sets | all UNSAT on both solvers, `--frames` |
| `cases_k6_expI.json` (23) | the last 23 k=6 cases, closed earlier by Lemmas B′/C | all UNSAT (slack ≤ 2 cliques, no `--frames`) |

- All 53 k=7 survivors have capacities (2,2,1⁵), and their two capacity-2 normals span an A2.
- (D)+(X) alone is UNSAT for 5 of the original 7 cases. (2,3,4,5,56,87,112) and
  (0,1,28,30,58,70,73) also need the frame clauses.
- The engine propagates zeros and never branches on supports. Branching over supports with SAT
  is what closes these cases.
- **Control.** On random exact integer decoders with true supports fixed (sub-universes of 40 or
  70 lines, including ± cancellation pairs and capacity-2 groups), the same clause families
  answer SAT in all 22 trials.

## Run

```sh
python3 check_S.py cases_k7.json --frames     # needs python-sat; ~2 min
python3 check_S.py cases_k9a.json --frames    # ~10 min each
python3 check_S.py cases_k6_expI.json
```

## Cost of a full k = 7 closure

This has now been done; see `../k7/README.md`. The estimate that used to be here (27.93M raw
(set, composition) pairs, about 3.2M after orbit dedup) came from an incomplete stage-1 run. The
complete run has 32,010,853 raw filter-surviving pairs and 3,636,308 after orbit dedup. All 7,536,480
(orbit, composition) cases were processed, which took about 11 CPU-h with the C clique scan. That left
2,342 engine survivors, which is more than the 53 sampled above. Not all of them are of the
(2,2,1⁵) A2 type: 222 have orthogonal capacity-2 normals and 9 have capacities (3,1⁶). With
`--frames`, 2,332 of them are UNSAT. The other 10 need wider (R) clique sets: slack ≤ 2 closes 8 of
them, and slack ≤ 4 closes the last 2.
