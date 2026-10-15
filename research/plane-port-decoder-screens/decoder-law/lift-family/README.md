# expSM: lift family of rank-9 masked completions on P120, cost ≥ 72 in any GL9 gauge (local; no κ claim; nothing posted)

Changed hypothesis (vs expSJ §4.2): stop treating the column side and the cross pairing numerically. Instead, write every
lift in one explicit symmetric normal form (the "κ-family") and bound the exact credit Σ(9 − ρ_k) of every dual-basis
decoder by a label argument. **Result: Theorem 1 below closes every lift with |S| ≥ 2, plus the concurrent A2 triples,
at cost ≥ 72 in arbitrary GL9 gauge.** Still open: |S| = 1 with special B, non-concurrent A2 triples, and non-lift
completions. **Stop test: no cost-≤71 configuration exists in any closed case. Every relaxation that reached credit 10
was traced to a modelling bug or a dropped true constraint (controls), so there was nothing to build in Fractions.**

## 0. Corrected parametrisation (checked against sols_0_20.npz: `chkparam.py`, `chkparam.log`)
Gauge: X_T = [R_T, 0] (torus-normalised; R = true root coordinates, r·r = 2). Then:
* x_s = ξ_s(c_s, 1);
* y_t = (r_t/2, η_t);
* y_s = (β_s r_s, η_s), i.e. the first 8 coordinates of y_s are β_s r_s, as expSJ had.

The constraint is η_t = −c_s·r_t/2 for t ⊥ s. The correction:
* (i) **c_s = c0 + μ_s r_s with a common c0** (concurrent) whenever the roots of S are independent. Pairwise,
  c_s − c_s' ∈ span(r_s, r_s'), because orth(s) ∩ orth(s') spans the complement. Independent directions then make the
  lines c_s + K r_s concurrent.
  * Numerically, the residual is ≤ 1e-13 on all 23 lift completions.
* (ii) Gauge c0 = 0. Then **η_t = 0 on T ∖ F, and η_t is free on F** := {lines non-orthogonal to all of S}.
  * Y_T has rank 8 iff η vanishes on F. |F| = 0, 1 or 0, 3, 6, 10, 16, 28, 56 for |S| = 8, 7, 6, 5, 4, 3, 2, 1
    (`cliques.py`).
  * The η_t on F are the extra freedom that expSJ's formula omitted.
* (iii) **The degenerate case ξ_f = 0 with a_f·r_f = 0 is vacuous.** If ξ_f = 0, then N_ft = 0 on orth(f) forces
  a_f ∝ r_f; then a_f·r_f = 0 gives x_f = 0, contradicting N_ff = 1. With T maximal, ξ_s ≠ 0 on all of S.
* (iv) **Exception: the A2 triples** (|S| = 3, r_s3 = ±r_s1 ± r_s2; 1120 sets, one W-orbit, F = ∅). Here c_s need not
  be concurrent: there is 1 extra parameter, and Y_T then has rank 9 (`a2triple.py`).

**κ-family (Prop. 0).** For S with independent roots, every lift is GL9 × torus-equivalent to:

    X = [diag(α) R, κ],   Y = [diag(α') R, κ']

with:
* supp κ = A = S and supp κ' = B ⊆ C(S) = S ∪ F;
* no a ∈ A, b ∈ B with a ⊥ b;
* 2α_aα'_a + κ_aκ'_a = 1 for every line a.

The family is symmetric in X ↔ Y, and α (resp. α') vanishes only on A ∩ B.

Census remark: 23 of the 28 saved completions are lifts (one side has rank 8 on T). **5 are not:** 0_1 and 13_1 are the
sign-twisted ones; 2_1, 5_1 and 7_1 have rank 9 on T on both sides. All 5 have numerical mincost 72 (expSJ).

## 1. Deficient-line classification (Lemma 2; exact control `kfam_ctrl.py` → PASS 25 + PASS 45)
Clean(K) means: for every class O, either rank_Q R[O∖K] = 8, or O = E7A1(f) with f ∈ K and kernel ∝ r_f.

If Clean(K) holds for a side with κ-support K, that side's deficient lines are exactly two kinds:
* **e9**, with Def = span(q_K)^⊥;
* **g_f = (κ_f r_f, −2α_f)**, with Def = ⟨q_f⟩, for each f ∈ K that is odd to all of K or has α_f = 0.

Proof:
* A direction (u, γ) vanishes on O ∖ K iff u ⊥ R[O∖K], so rank 8 forces u = 0.
* The γ = 0 case has the same argument.
* The Def of the standard lines uses min_{O≠E7A1(f)} |O ∖ E7A1(f)| = 24 (`p1ext.py`).

Finite facts:
* Clean(K) holds for every |K| ≤ 13, since |O| ∈ {56, 64}, hyperplanes other than E7 hold ≤ 42 lines, and m1 = 24.
* **Clean(C(S)) holds on every clique W-orbit with |S| ≥ 3** (`orbits.py`, `orbits.log`):
  * the orbits number 1, 1, 2, 1, 1, 1, 2, 1 for k = 1..8;
  * |C(S)| = 57, 30, 3/19, 14, 11, 9, 7/8, 8.
* Clean is inherited by subsets B ⊆ C(S).
* For |S| = 2 exactly 3 classes are non-Clean (h = 499, 500, 504; each contains all of F).
* For |S| = 1 all 509 are non-Clean.

## 2. Theorem 1 (cross pairing; hand proof, label-only)
**Statement.** For any κ-family member with A, B ≠ ∅ and Clean(A), Clean(B), every square decoder (any dual-basis pair,
any GL9 gauge) has Σ_k (9 − ρ_k) ≤ 9, so **cost ≥ 72**.

Here ρ_k is exact: ρ = (9 − dim D_R) − dim(D_R^⊥ ∩ D_A).

Pairings:
* e9·e9' = 1;
* e9·h'_c = −2α'_c;
* g_f·e9' = −2α_f;
* **g_f·h'_f = 2(κκ' + 2αα') = 2 ≠ 0**.

H_X ∪ H_Y (the indices of the standard lines) is an orthonormal label set: A × B pairs are odd, so q·q = 0. Hence it is
independent and has ≤ 9 elements.

Cases:
* **Neither e9 nor e9' used.** Each position earns ≤ (number of distinct indices it uses), and g_f, h'_f must share a
  position. Total ≤ |H_X ∪ H_Y| ≤ 9.
* **Both used.** e9 and e9' are forced into one position, with credit 9 − d_A + dim(U_A ∩ U_B^⊥).
  * The other positions may only hold g_f / h'_c with α = 0, so with indices in A ∩ B. Their labels W2 ⊆ U_A satisfy
    W2 ∩ U_B^⊥ = 0.
  * Total ≤ 9 − d_A + d_A.
* **Only e9 used** (only e9' is symmetric). The Y-line at e9's position is some h'_c, with credit 9 − d_A + [q_c ∈ U_A].
  * g_c is excluded.
  * The other indices J ⊆ H_X ∖ {c} are orthonormal in U_A and orthogonal to q_c (q_c·q_c = 1).
  * Total ≤ 9.

∎

**Interpretation of the requested statement.** Terms deficient on both sides do exist, e.g. (g_f, h'_f) with credit 1,
or (g_f, h'_f') with credit 2. They never gain, because the forced non-zero pairings (the diagonal identity
2αα' + κκ' = 1 and e9·e9' = 1) pin their partners.

**Corollary (proved).** Assume a lift with T = V∖S rigid of rank 8 and S pairwise odd. Then C ≥ 72 in arbitrary gauge
in each of these cases:
* **(a)** every |S| ≥ 3 that is not an A2 triple (all B, all special parameters). This includes |S| = 8.
* **(b)** concurrent A2 triples;
* **(c)** |S| ∈ {1, 2} whenever B is Clean;
* **(d)** |S| = 2 with non-Clean B (B ⊇ F minus ≤ 1 line). `relax2.py` handles this case:
  * Extra Y-lines all lie in W3 = P × K (P = A2 plane): at most 3 of them, with γ' ≠ 0.
  * Their directions are among the 6 special P-directions or are generic.
  * Def is upper-bounded (credit is monotone in Def).
  * Lines carrying h499 form a subspace of dimension ≤ 2.
  * Zero patterns are chosen by CP-SAT. Result: **max credit 9 on all 8 non-Clean B, PASS.**

So for |S| = 2, every B is closed.

## 3. Checks (all `../expO/run <script>`; logs in this dir)
| script | content | result |
|---|---|---|
| `chkparam.py` | §0 normal form on sols_0_20 | all lifts fit; 5 non-lifts flagged |
| `kfam.py`, `kfam_ctrl.py 11 25` | exact random κ-members incl. α/α'=0 patterns, \|S\| = 1..8: exact completion, exact deficient lines = Lemma 2 prediction, pairing identities | PASS 25 (and 45 earlier) |
| `p1ext.py` | m1 = m2 = 24; \|O∩orth(s)\| distribution | — |
| `orbits.py` | clique W-orbits; Clean(C(S)) per orbit | as §1 |
| `relax.py`, `run_s8.py 8` / `7` | CP-SAT exact credit with forced/free pairings, all B ⊆ S and zero patterns (S_k-symmetric, isometry argument) | max credit 9 (495 / 330 cases) |
| `run_orbits.py 1 60` | the engine on random Clean B for every orbit k = 2..8 | PASS 480, max 9 |
| `relax_ctrl.py` | **negative controls**, each dropping one true constraint: an infeasible α = 0 outside B; g_f·h'_f may vanish; e9·e9' may vanish | 10 / 16 / 10, and the true model gives 9: **PASS** |
| `relax2.py` | \|S\| = 2 non-Clean B | PASS, max 9 |
| `a2triple.py 3 4` | A2 triples, exact: concurrent and non-concurrent | mincost 73 (generic) |
| `k2special.py 5 2` | stop-test probe: \|S\|=2 special locus (E1 exists; E1 merged with h'_s1, Def 2-dim) | mincost 72 |

Caveats:
* `relax2.py` passed only after two fixes. The fixes were a name-collision bug, and twin variants for Def bounds
  containing h499; an earlier looser version gave 10–12 for exactly those reasons.
* The hand proof (§2) is the primary argument for (a)–(c); the engines are confirmations.

## 4. What remains
1. **|S| = 1 with non-Clean B.**
   * Here B ⊆ C(s) has 57 lines. The extra Y-lines form the pencil (λ r_s, γ'), with zero sets equal to the level sets
     of κ'_b/(α'_b r_b·r_s) on nonorth(s), plus lines with u' ⊥ the rank-≤6 sets O ∩ orth(s).
   * Neither Theorem 1 nor `relax2.py` covers this case.
   * Next: the same W3-type engine, with the 1-dim pencil and with Def sets indexed by level-set partitions. Note that
     h'_s, e9 and the pencil span only 2 dimensions.
2. **Non-concurrent A2 triples.** These have 1 extra parameter, and Y_T has rank 9. The generic exact mincost is 73;
   special parameters are not analysed.
3. **Classification of completions.**
   * Rigidity of T, and the claim that every rank-9 completion is a lift, are not proved.
   * 5 of 28 saved completions are not lifts in the sense used here: the sign-twisted ones, plus 2_1, 5_1 and 7_1.
   * A tangent-dimension comparison at one point per lift type (expSJ next-experiment 2) is still the cheapest test.
4. m ≥ 10 terms and non-square decoders: untouched.

## 5. Next experiment (ranked)
1. **|S| = 1 closure.** Extend `relax2.py` to K = C(s): enumerate level-set partitions of nonorth(s) up to Stab(s),
   using Def upper bounds.
   * *Stop test:* any credit-10 configuration that survives the exact pairing check is built as an exact Fractions
     member and decoder, then checked.
2. **Tangent dimension at points 2_1 / 5_1 / 7_1 and 0_1** vs. the lift-family dimension, to decide whether they lie
   on new components.
   * *Stop test:* a new component with a numerical mincost below 72.
3. Symbolic rank jumps for the A2 non-concurrent family. It is a 1-parameter family, so sympy is feasible.
