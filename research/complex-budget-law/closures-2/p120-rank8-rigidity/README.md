# P120 rank-8 rigidity and sparse-transfer ports: idea 4 and the corrector-rank lever closed for line ports at h ≤ 12

Status: scoped closure. The rigidity theorem has a machine-checked finite step, and the sparse-support screen has exact lower
bounds. Author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`. Multipliers are over κ = 0.0017.
This lane merges two scouting runs: expNO (floors and half-integral decoders) and expQF (rigidity and the exact sparse
screen).

## Verdict

**Multiplier: ×1.00.** Two levers are closed for line ports at h ≤ 12:
- idea 4, sparse-transfer ports;
- the expNT corrector-rank lever. It stands or falls with idea 4, because r(a) ≤ d_T(a).

| Result | Status |
|---|---|
| **Rigidity theorem.** Over any field of characteristic ≠ 2, every masked decoder of rank ≤ 8 on the E8 P120 labels equals Λ(G/2)Λ⁻¹ for an invertible diagonal Λ. Hence T = F (all 56 free pairs per port are transfers), N_ab N_ba = 1/4 on free pairs, and max r(a) = 8 = h − 1. | proved; Step 1 is a finite Boolean check, certified by UNSAT from three solvers plus our own DRUP/RUP check |
| **Exact sparse screen.** Every structured transfer pattern from E8 sub-root systems with d̄ ≤ 30 has certified min-rank ≥ 12. That covers the singles A8, D8, E7+A1, A4² and D4², and 27 unions of Weyl conjugates. Every union with d̄ ≤ 30 has ≥ 15. | exact lower bounds (identity principal minor on an independent set of T) |
| **Half-integral decoders.** T = F exactly, for definite or indefinite W. | proved; checked on E8 and Z^{9,1} |
| **Universal transfer floor.** d̄_T ≥ about 1/γ − 1 (0.75 at γ = 4/7); with full-cost terms, h/γ − 1 (14.75). | proved (LP), so counting alone does not force density |
| **Caro–Wei-tight configurations.** γ ≥ 1. | proved |

## Changed hypothesis tested

E8 routes every free pair as a lossless transfer, and helpers pay for transfers. The hypothesis is that this density is
incidental: some compressing decoder of rank ≤ 10 (line ports, h ≤ 12) could use a sparse T (d ≤ 26–30). If so, the bound
r(a) ≤ d_T(a) would open the corrector-rank lever as well.

## Rigidity theorem: proof

**Setup.** N = XYᵀ with X, Y ∈ K^{120×8}, char K ≠ 2, N_aa = 1, and N_ab = 0 for every orthogonal pair a ⟂ b.

**Step 0 (normal form).**
- Fix an orthogonal frame Q = {f_1, …, f_8}. Then N[Q,Q] = I.
- Use the GL₈ gauge to set X_{f_j} = e_j and Y_{f_j} = e_j*.
- Every other root is a = ½ Σ_{j∈S_a} ε_j f_j, where S_a is a block of the Steiner system S(3,4,8). Each of the 14 blocks
  carries all 16 sign patterns.
- Since a ⟂ f_k for k ∉ S_a, both X_a and Y_a are supported on S_a.
- For blocks S ≠ S′ with S ∩ S′ = {i,j}:
  - a ⟂ b ⇔ the types t_ij(a) = ε_iε_j and t_ij(b) differ;
  - N_ab = X_a(i)Y_b(i) + X_a(j)Y_b(j).

**Step 1 (no zeros; finite check).**
- Keep only zero patterns, with two conditions: an orthogonal pair's sum cannot have exactly one nonzero term, and N_aa = 1
  needs at least one nonzero term.
- These conditions together with "some entry X_a(j) or Y_a(j), j ∈ S_a, is zero" are **UNSAT**.
- Certificates:
  - `nondeg_cnf.py`: 13,442 variables, 49,874 clauses; CaDiCaL and Glucose both return UNSAT;
  - `rupcheck.py`: our own RUP checker verifies the 2,144-lemma Glucose DRUP refutation;
  - the scouting run also checked it with CP-SAT, both with and without the 2,025 frame-invertibility constraints.
- A deletion-minimal core still uses 31 roots, so no short hand proof of this step is known.

**Step 2 (each pair is rigid in 2D).**
- Fix a pair {i,j}. Write r_a = X_a(j)/X_a(i) and s_b = Y_b(j)/Y_b(i); these exist and are nonzero by Step 1.
- For orthogonal a (type τ) and b (type −τ) in different blocks, r_a s_b = −1.
- The cross-block bipartite graph is connected (K_{3,3} minus a matching, blown up).
- So r_a = λ^τ_ij for every root of type τ, and s_b = −1/λ^{−τ}_ij for every root of type τ as well.

**Step 3 (frame inside a block).**
- In a block, ½(++++) ⟂ ½(++−−).
- N_{aa′} = 0 then reads −1 − λ⁺_ij/λ⁻_ij = 0, so λ⁻_ij = −λ⁺_ij =: −w_ij.

**Step 4 (Steiner cocycle).**
- Reading X_{++++} along the triple ijk gives w_ik = w_ij w_jk.
- Every triple lies in exactly one block, so w_ij = d_j/d_i.

**Step 5 (conclusion).**
- X_a = α_a D c_a and Y_a = β_a D⁻¹ c_a, with c_a(j) = (a·f_j)/2.
- So N_ab = α_aβ_b (a·b)/2, with α_aβ_a = 1. ∎

**Gram reading.** Up to projective equivalence, the E8 roots are the only realisation of their own orthogonality relation as
a point/hyperplane incidence in P⁷.

**Local cross-check (`tangent.py`).** The Jacobian of the decoder equations at G/2 (7,680 × 1,920) has exact rank
1,737 = 1,920 − 183, and the 183 gauge directions lie in its kernel. So the variety is smooth of dimension 183 at G/2 and
equals the gauge orbit there.

## Exact sparse-support screen

**Pattern.** For a sub-root system R ⊂ E8, T_R = {(a,b) ∈ F : a ∓ b ∈ R}.

**Certified bounds.**
- Lower bound: α(T) per component. N restricted to an independent set of T is the identity, so rank N ≥ α(T).
- Upper bound: block G/2, a rational completion, where T = F inside the block. Otherwise there is no rational upper bound,
  and `alpha.py` prints a sentinel above 999.

| Pattern | d̄ | min-rank |
|---|---|---|
| D8 | 26.1 | = 16 (exact) |
| A8 | 16.8 | ≥ 16 |
| E7+A1 | 29.9 | ∈ [12, 16] |
| A4² | 9.3 | ≥ 24 |
| D4² | 11.2 | = 32 (exact) |
| 27 unions of Weyl conjugates with d̄ ≤ 30 (`unions.py`, seed 0, 3 reps) | 17–30 | ≥ 15 (A4²∪A4²: 23; D4²∪D4²: 32; any union with A8 or D8: 15–16) |

**Heuristic pricing (`price_rank.py`).** This uses expNO's extrapolation γ ≈ 4/7 + 0.075(K − 8).
- K = 12, d = 30: γ ≈ 0.87, which prices at pair ×0.06 and Λ ×0.08.
- K ≥ 15: γ > 1, so there is no compression.

## Floors and the half-integral lemma (from expNO)

**(a) Universal floor (`floor.py`).**
- Caro–Wei on T gives α(T) ≤ n₊, the number of positive-cost totals, so d̄ ≥ V/n₊ − 1.
- Hadamard capacity gives V ≤ Σ(2^{β_j} − 1), and C = γV.
- The LP is solved exactly by vertex enumeration: d̄ ≥ about 1/γ − 1.
- With generic full-cost terms the floor is h/γ − 1. E8's 52–56 sits about 3.5× above it.

**(b) Caro–Wei-tight configurations never compress.** Reader columns are T-cliques, so C ≥ V.

**(c) Half-integral decoders (`halfint.py`).**
- Setup: W is integral with W_aa = 2 and off-diagonal entries in {0} ∪ odd.
- Then N = W/2 is a valid decoder for B = J + (W mod 2), and h ≤ rank W + 1.
- A free pair has W_ab odd, so N_ab ≠ 0. **Hence T = F.**

**Screens (float, not included).** E8 subsets with d ≤ 26 give γ ≥ 1.5. Indefinite sets in Z^{9,1} give γ ≥ 1.09 at
d ≤ 24, and the unconstrained optimum is γ = 0.849 at d̄ = 50. Rank-9 and rank-10 sparsification of the E8 labels stays dense
(d̄ = 55.5 and 49.5).

**Payoff table (`score.py`).** Even a hypothetical sparse decoder at E8's γ is worth only ×1.0–1.9 (Λ model) or ×1.01 (pair
model).

## Stop test that fired

The test was: stop if min-rank ≥ 12 for every structured pattern with d ≤ 30. It fired:
- singles ≥ 12;
- unions ≥ 15;
- rank 8 is rigid.

Neither continue test (rank ≤ 10 with d ≤ 26, or rank ≤ 8 with d ≤ 40) can fire.

## Not closed (each would need a new idea)

- Supports at rank 9–11 with 30 < d < 50 that are not root-structured. Only LM probed these, and the fits stayed dense.
- Labels that are not E8.
- **Reopen test:** a rank-≤10 masked decoder with d̄ ≤ 26, γ ≤ 0.6 and a certified α(T) ≤ 10.

## Files and reproduction

Run from this directory. Runtimes are for one worker under a 4G cap.

| Script | Run | Runtime | Output |
|---|---|---|---|
| `rupcheck.py`: RUP check of `nondeg.drup.gz` against `nondeg.cnf.gz` (stdlib) | `python3 -I -B rupcheck.py` | 23 s | `rupcheck.out`, `PASS_RUP` |
| `nondeg_cnf.py`: rebuilds the Step 1 CNF, checks it is identical to the stored file, and runs CaDiCaL and Glucose (`--write` regenerates the CNF and proof) | `uv run -q --no-project --python 3.12 --with numpy --with python-sat python -I -B nondeg_cnf.py` | 3 s | `nondeg_cnf.out`, `PASS_NONDEG_CNF` |
| `tangent.py`: exact tangent-space rank. **Needs python-flint.** | `uv run -q --no-project --python 3.12 --with numpy --with python-flint python -I -B tangent.py` | 24 s | `tangent.out`, `PASS_TANGENT` |
| `alpha.py`: exact α per component for the singles (CP-SAT) | `uv run … --with numpy --with ortools python -I -B alpha.py` | 2 s | `alpha.out`, `PASS_ALPHA` |
| `unions.py`: union screen | `uv run … --with numpy --with ortools python -I -B unions.py 0 3` | 42 s | `unions.out`, `unions_0.json`, `PASS_UNIONS` |
| `floor.py`: exact floor LP (stdlib) | `python3 -I -B floor.py` | 0.1 s | `floor.out`, `PASS_FLOOR` |
| `halfint.py`: half-integral T = F (stdlib) | `python3 -I -B halfint.py` | 0.7 s | `halfint.out`, `PASS_HALFINT` |
| `score.py`: payoff table (stdlib), using `rho_min.py`, `kappa_pair.py` and `transfer.py`, copied with only their import paths changed | `python3 -I -B score.py` | 0.2 s | `score.out`, `PASS_SCORE` |
| `price_rank.py`: heuristic pricing of the certified patterns | `python3 -I -B price_rank.py` | 0.1 s | `price_rank.out`, `PASS_PRICE` |
| `e8.py`, `subsys.py`: the P120 model and sub-root systems | — | — | — |

Not included from the scouting runs:
- `nondeg.py` (the CP-SAT form of Step 1);
- `core.py`;
- `globnum.py` (5 of 8 LM starts converged, all to the G/2 orbit);
- `screen1.py`;
- the expNO float screens.

Prepared with Claude Code (Claude Opus 5.5) assistance. Not independently reviewed.
