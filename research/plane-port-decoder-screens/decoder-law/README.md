# expSG: the "C = |L| − 48" law for masked line decoders (local; no κ claim, nothing posted)

Model: L ⊆ the 120 E8 root lines. Decoder N = Σ_k v_k w_kᵀ over Q with N_aa = 1, N_ab = 0 when a ⊥ b,
free when ip = ±1. Term cost ρ_k = rank_F2 B[supp v_k, supp w_k], with B = (J + G) mod 2. This equals the
label pairing q_a·q_b; `e8.py` asserts it on all 14,400 pairs. C(L) = min Σρ_k and gain(L) = |L| − C(L).

## 1. Main result: the law is one statement about P120 (exact, elementary)

* **Restriction.** If L' ⊆ L then C(L') ≤ C(L). Restrict every term; ranks of submatrices can only drop.
* **Blocks and padding.** C(L₁ ⊔ L₂) ≤ C(L₁) + C(L₂) for **any** disjoint L₁, L₂ (no orthogonality needed),
  because block-diagonal N is always valid: a cross entry 0 is allowed on both required and free pairs. Also
  C({a}) = 1 (B_aa = 1).
* Hence gain is **monotone and superadditive**: gain(L') ≤ gain(L) for L' ⊆ L, and
  gain(L₁ ⊔ L₂) ≥ gain(L₁) + gain(L₂).
* **Corollary.** "C(L) ≥ |L| − 48 for every E8 line set L" ⇔ **C(P120) ≥ 72**.
  * Any counterexample anywhere (P112 ≤ 63, P111 ≤ 62, P104 ≤ 55, any sub-configuration with gain ≥ 49) pads
    with identity terms to a P120 decoder of cost ≤ 71.
  * Conversely, C(P120) ≥ 72 proves P112 ≥ 64, P111 ≥ 63 and P104 ≥ 56 together.
* The **complex** lever needs C(P112) ≤ 56, so it needs C(P120) ≤ 64. The **bit** lever needs C(P111) ≤ 56, so it
  needs C(P120) ≤ 65. These are statements about the same single object. (This holds insofar as the bit decoder is
  the same masked line problem on P111. That is the scout's framing; it is not re-derived here.)
* The upper bounds are known: P120 frame decoder 72, P112 64, P111 bit 63. So gain = 48 exactly on every L
  between P112 and P120 **iff** the law holds.
* `verify.py` (stdlib) checks all of this exactly over Fractions:
  * the P120 frame decoder is exact (72);
  * its restriction to P112 is exact (64);
  * padding back gives 72;
  * random restrictions and paddings, and a two-frame block decoder, are exact.
  It prints PASS.

**Status of C(P120) ≥ 72 = the law: open in general.** It is at least as hard as the P112 64 floor.

## 2. Proved / exhaustively checked restricted classes on P120

(a) **Decoders with ≤ 8 terms cost exactly 72.** Proved, modulo 8-term rigidity, which is checked numerically.
* rank N ≥ 8, because a frame gives N[F,F] = I₈.
* A rank-8 completion on P120 is a torus orbit D(G/2)D⁻¹. LM finds this: 6/6 sym/asym runs, residual
  < 2e-14, |N_ab N_ba| = 1/4 on all odd pairs (`rig8.py`). It is the P120 analogue of #1212.
* Given that, every term is v = D R m, w = D⁻¹ R m'. Its support contains Φ ∖ H for a rank-7 flat H.
* All 8 ray types (W-orbits of fundamental weights) have supports whose labels span F2⁹ (`verify.py`).
* The form on F2⁹ is nondegenerate, so ρ = 9 for every term. Hence C = 72.
* This also covers every "lattice gauge" (any basis/dual basis M, M⁻ᵀ) decoder on P120.

(b) **Exhaustive Boolean-support SAT** (`exsat.py`; CaDiCaL via pysat, no proof logs). See the table in §4.

## 3. Searches for gain > 48 (constructions are exact when found)

* Lattice-gauge decoders (G/2 in an arbitrary ray basis), annealed jointly over the basis and the subset L
  (`gsearch.py`, 2 × 40k steps): max gain **48**.
  * On P120 every basis gives exactly 72, as (a) predicts.
* Block decoders, 2–3 blocks with independent gauges plus identity (`part.py`): max gain **48**. The best state is a
  112-line block with gain 48 plus one singleton.
* Sub-root-systems, best gauge (`gfix.py`): gain is 0 or negative on every proper subsystem of rank ≤ 7 and on D8,
  A8, E6+A2, A4+A4, …. It is 14 on E7 and E7+A1, and 48 only on P112/P120.
  * Gain is concentrated, not ~6·rank. The structural form is "48 = gain(P120) = 120 − 8·9".

## 4. Exhaustive class results on P120 (Boolean relaxation; necessary for an exact decoder over any field)

| class (P120, budget ≤ 71) | result |
|---|---|
| control: 8 incumbent frame supports | budget 72 SAT, 71 UNSAT |
| ρ≤7 symmetric hyperplane sections + full supports | SAT (7 full + one ρ7) — relaxation blind to full supports; excluded algebraically by §2(a) (8 terms) |
| **only ρ7 symmetric hyperplane sections (255 classes, 1 copy each)** | **SAT at cost 70**: 10 terms h = 24,30,34,36,40,46,48,54,58,60 (5 E7+A1-type, 5 D8-type; h span a 4-dim W) |
| that certificate, LM 16 restarts sym/asym | **no exact decoder**: residual² = 51 (sym) / 50 (asym), integer plateaus as in the scout's P112 runs |
| random families (p120lm.py): 9×E7+A1, 10×E7+A1, 10×ρ7 mixed, 9 sym hyperplanes ≤71, 2,995 asym hyperplane pairs ≤71 | all Boolean-UNSAT (≈3,600 families) |
| only E7+A1 class / only D8 class, budget 71 | timed out (1400 s), inconclusive |

So the Boolean relaxation does not close the hyperplane class; the one surviving certificate fails algebraically
(numerically). The obstruction in [56,72) is algebraic, as the scout found on P112.

## Overlap (coordinator alert, read PR bodies only)

* **#2338 (DreamingOfClouds, open)** already proves E8 (=P120) ≥ 72, P112 ≥ 64, P104(#968) ≥ 56 for factors in
  col M of a rank-8 realisation (any m), and "≤ 8 paid terms ⇒ ≥ 64"; Theorem R (dense B at h = 9). My §2(a) is
  a re-derivation of the E8-72 / ≤8-term case (numerical rigidity route) — credit #2338.
* **#2724 (DreamingOfClouds, open)** closes m = 9 case (b) (one factor outside col M, rank-8 realisation) on
  P112 (≥ 64) and P104 (≥ 56). Open there: case (c) square decoders of rank-9 realisations (co-star families),
  m ≥ 10, h ≥ 10, other punctures.
* **#2719 (ywz65)**: D_n root-code all-term floor n(n−1); a sub-system floor, consistent with gain ≤ 0 on D_n here.
* Not seen in those bodies: the monotone/superadditive gain reduction (§1). With it, #2338's E8-72 and its
  P112-64 are the same statement once extended beyond col M; any future general E8 ≥ 72 proof yields P112 ≥ 64,
  P111 ≥ 63 and every puncture floor |L| − 48 for free, and the uncovered cases (rank-9 square, m ≥ 10, h ≥ 10)
  can be attacked on P120 alone, where W(E8) acts transitively and every term has ρ ≤ 9.
* My §4 certificate is an m = 10, h = 9 (dense B) family — inside #2338/#2724's OPEN "m ≥ 10" region; it fails
  numerically.

## 5. Implications

* Decoder/centre lever (complex P112 and bit P111): **closed iff C(P120) ≥ 72**. Every route to a cheaper
  decoder on either side must produce a P120 decoder at ≤ 71 (≤ 64 for the ×1.22 complex target).
* Rigidity + (a) mean that such a decoder needs ≥ 9 terms with average ρ < 8. Therefore some terms have supports
  whose labels lie in a hyperplane: they omit ≥ 56 of the 120 lines on one side.

## Next experiment (ranked)

1. Exhaustive Boolean SAT on P120 for all symmetric and asymmetric hyperplane / codim-2 section pairs, with
   multiplicity, at budget 71.
   * If UNSAT, the law holds for every decoder whose terms all avoid full-span supports.
   * The remaining class then has ≥ 1 full-span (ρ = 9) term, which the Boolean relaxation cannot see. That class
     needs the rank-9 single-lift analysis.
   * Stop test: a SAT certificate followed by LM failure (≥ 8 tries per family, ~25% control hit rate) on ≥ 50
     certificates.
2. Prove 8-term rigidity on P120 exactly (deficiency-space method of PR1248/PR1497). This turns (a) into a theorem.

Note: LM on P120 recovers the incumbent in only ~1 of 4 random starts. Single-start residuals around √50
are not evidence of infeasibility, so the scout's P112 residual statistics should be read with ≥ 8 tries.
