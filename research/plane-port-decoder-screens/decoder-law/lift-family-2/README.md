# expL1: closing the rank-9 square lift family for C(P120) ≥ 72 (local; no κ claim; nothing posted)

Builds on expSM (κ-family normal form, Lemma 2 = Clean classification, Theorem 1 = credit ≤ 9), expSJ (census,
`maxdef.rho_def` / `mincost`) and expSG (reduction: every decoder floor ⇔ C(P120) ≥ 72). Credit of a position is
9 − ρ = dim D_X + dim(D_X^⊥ ∩ D_Y); cost ≤ 71 needs total credit ≥ 10 over the ≤ 9 positions of a dual-basis pair.

**Changed hypothesis.** The open cases of expSM §4 do not need a new search. Each has a finite *line geometry*:
* |S| = 1: an E7 frame plus the 56 weights, with 27-sets M_b^±.
* A2 triples: three 2-dim kernels with pairwise intersections.

Once that geometry is computed exactly, the credit bound is a small exact optimisation over labels, which runs in
minutes. Census correction: the 5 "non-lift" completions of expSJ/expSM are lifts too, with hidden |S'| = 7 or 8.

**Result.**
* Every κ-family lift is closed at credit ≤ 9 (cost ≥ 72 in arbitrary GL9 gauge). This covers every A, B and every
  special parameter, including the cases left open by expSM:
  * |S| = 1 with non-Clean B;
  * non-concurrent A2 triples, all parameters;
  * a gap not listed in expSM: non-clique independent supports.
* All 28 saved completions are κ-lifts, and at each one the tangent dimension equals the κ-family dimension.

**What remains** is stated in §5. **Stop test:** no credit-10 configuration survived anywhere. Every engine has a
negative control that reaches ≥ 10 when one true constraint is dropped. Every exact member tested (86) has an
exact gauge floor of 72 or 73.

## 1. |S| = 1, every B ⊆ C(s), every level structure — CLOSED

Setup (s = line 0; E = orth(s), 63 lines; N = nonorth(s), 56 lines; C(s) = {s} ∪ N):
* Fix the sign so that r_b·r_s = 1, and write r_b = r_s/2 + p_b. The p_b are the 56 weights of E7, with
  2p_b·p_c ∈ {3, 1, −1, −3}.
* The level of b is L_b = κ'_b/(α'_b r_b·r_s), and level(s) = κ'_s/(2α'_s), which is ∞ if α'_s = 0.
* Write a Y-line as F = (u0 ∈ r_s^⊥, t, γ). F vanishes:
  * at a ∈ E iff r_a·u0 = 0;
  * at b ∈ N iff p_b·u0 + t + γL_b = 0;
  * at s iff 2t + γL_s = 0.

**Lemma A (Y-line classification; exact, `geom1.py`, `geom2.py`, `facts.py`, control `kctrl1.py`).** Every Y-deficient
line is one of two kinds.
* **Pencil line** (u0 = 0). It vanishes on E ∪ (one level block P), with **Def = {h : O_h ∩ C(s) ⊆ P}**.
* **Generic line**, with Def ⊆ {H0, h_b, h_b*}. These classes are:
  * H0 = 504, with O = N;
  * h_b, with O_{h_b} = (E ∩ p_b^⊥) ∪ {b} ∪ M_b^−, where M_b^− = {c : 2p_b·p_c = −1} (27 lines), and h_b ⊕ h_b* = H0.

  When these classes occur:
  * **h_b occurs iff L is constant on M_b^− and L_b differs.** The line is then unique.
  * **H0 occurs iff L is affine and non-constant on N.**

Proof. The only classes with rank R[O∩E] ≤ 6 are H0 and the 56 classes h_b. Their E-parts force u0 ∝ p_b (`geom1`).
Solving the remaining conditions gives the stated criteria; γ = 0 is impossible for a generic line. Kernels are always
≤ 1-dimensional when rank Y = 9. ∎

**Finite facts** (`facts.py`):
* Non-antipodal M's always intersect (10 or 16 common lines). Therefore all constant M's lie in one block, or in two
  antipodal blocks.
* The largest set of weights with pairwise 2p·p = 1 has 7 elements.
* The labels of C(s) span F2⁹.

**Case I: X-line e9 used (hand proof).**
* Every other position has a Y-line ⊥ e9, so γ = 0. Such a line is deficient only if it is (r_s, 0) with α'_s = 0,
  and then Def = ⟨q_s⟩.
* The partner of e9 earns 8 + [q_s ∈ D_Y]. This holds only for h'_s, which needs α'_s ≠ 0.
* If g_s is also used, every other line is (u0, 0, 0), which is non-deficient. The partner of e9 cannot be h'_s,
  because g_s·h'_s = 2 ≠ 0.
* So **credit ≤ 9**.

**Case II: neither e9 nor g_s used.** Credit = T_Y, the total Y-deficiency of the used lines.
* `s1relax.py` is a CP-SAT label model. It has:
  * pencil blocks 1–2 (at most 2 pencil lines, since they span a 2-dim space);
  * ≤ 2 further constant-M blocks;
  * a generic line per b with M_b^− in a block not containing b;
  * a cap of ≤ 3 lines in each W_b = {(λp_b, t, γ)};
  * Y rank 9;
  * Def dims encoded as |D| ≥ 2^j − 1.
* **Result: max T_Y = 9**, e.g. one 56-line block plus F_77. That configuration is the gauge image of B = {77}.
* The H0 case (`hz.py`): the pencil blocks are fibres of w, over 91 candidate directions plus all p_b. **T_Y ≤ 3.**

**Case III: g_s used, e9 not.** Credit = 1 + Σ dim D_Y − [the partner g'_0 is a pencil line whose block contains s].
* Every other Y-line lies in H* = {tκ_s = α_sγ}. Consequences:
  * The only pencil line in H* is ℓ*. Its level is λ* = −α_s/κ_s, and ℓ* never vanishes at s because 2αα' + κκ' = 1.
  * A generic F_b lies in H* iff λ* = (3v + L_b)/4, where v is the level of M_b^−. So such b with the same M-block
    share one level, which is neither v nor λ*.
* `s1gs.py`: **max credit = 9** (OPTIMAL). The binding configuration is a 51-line block containing s, with Def dim 3,
  plus 6 generic h_b lines: 1 + 3 − 1 + 6 = 9.

**Controls.**
* Negative controls, each dropping one true constraint, all reach ≥ 10:
  * `s1gs.py noS` (drop the q_s correction of D_0) → 10;
  * `s1gs.py sStar` (allow s in the ℓ* block) → 10;
  * `s1relax.py nob` (allow b in its own M-block) → 16.
* `s1gs.py ctrl` (drop the H*-level constraints) stays at 9, so those constraints are not binding.
* `kctrl1.py` builds exact Fractions members of 6 structured kinds: generic levels, Q-configurations with up to 7
  lines, antipodal blocks, affine along p_b, affine along an E7 root, and the engine's extremal configuration.
  * For each, it computes **all** deficient lines on both sides by exact kernels over the 510 classes and checks
    Lemma A.
  * **PASS 48/48.**
  * The exact-ρ gauge floor `mincost` = **72 on all 48**.

**⇒ Theorem L1. Every |S| = 1 lift (A = {s}, any B ⊆ C(s), any κ', any α'_s, including α'_s = 0) has credit ≤ 9.**
By X↔Y symmetry this also covers B = {b} with any A ⊆ C(b).

## 2. Non-concurrent A2 triples, all parameters — CLOSED

Normal form (`a2lib.py`): gauge c0 = 0, ν1 = ν2 = 0, ν3 = 1, using the GL9 shift and the 9th-column scale. Then:
* y_t = (r_t/2, −(n·r_t)/2) on T_3 (lines orthogonal to s3 only), and y_t = (r_t/2, 0) elsewhere on T, where n ∈ P ⊥ r_s3.
* **Y_T has no free parameter.** The only freedom is in μ_i, ξ_s and the directions of y_s.

**Lemma B (exact, `a2scan.py`, `a2scan2.py`).** Only six classes have a non-zero kernel on Y[O∩T]:
* E7A1(s_i), with 2-dim kernel K_i;
* the pair classes h_ij = q_i ⊕ q_j, with 1-dim kernel L_ij.

Moreover K_i ∩ K_j = ⟨L_ij⟩.

Hence the Y-deficient lines are exactly:
* h'_i, the unique line of K_i that vanishes at s_i. It is unique because y_si ∈ A_i ≠ 0.
* With at most one merge: **h'_i = h'_j = L_ij**, which happens iff y_si, y_sj ⊥ L_ij. Its Def is span(q_i, q_j).

The X side is e9 (Def dim 6) plus g_1, g_2, g_3 (⟨q_i⟩). g_3 ∈ span(e9, g_1, g_2), so at most 3 of these 4 can be used.

**Forced pairings.**
* **g_i·h'_i ≠ 0.** A_i × K_i is a perfect pairing, because RX_i × RY_i carries G/2 of rank 7. Also x_si ∈ K_i,
  y_si ∈ A_i and x_si·y_si = 1. If h'_i annihilated g_i, it would annihilate g_i and y_si, which span A_i.
* γ(L_ij) = 1, so **e9·L_ij ≠ 0**.

**Credit.** `a2credit.py` enumerates every position assignment exactly.
* **Max credit = 9**, PASS. Example: e9 alone (6) plus h'_1, h'_2, h'_3 alone.
* Negative control: dropping g_i·h'_i ≠ 0 gives 11.
* Exact members (`a2ctrl.py`, 32): generic, β_s = 0, merged pairs, and μ = 0 variants. Deficient lines match Lemma B,
  and **mincost = 73 on all 32**.

## 3. The 5 "non-lift" completions are lifts; tangent dimensions (`census3.py`, `tangent.py`, `hidden*.py`)

The expSJ/expSM "S" was defined by |N_ab N_ba| − 1/4 > 1e-6. That test misses lines of A∖B, and lines of A∩B with
tiny κκ'. A line s with κ_sκ'_s = 0 has **no** deviation against T, because N_sb N_bs = α_sα'_s/2.

Instead, take the support of the sparsest X-deficient line, then fit the κ-family with the c0 shift. Results:
* 0_1: |A| = |B| = 7. The κ are tiny (~1e-3); σ9(N) = 3e-6.
* 13_1: |A| = |B| = 8, also with tiny κ. **These two are not sign twists. They are near-rank-8 κ-lifts.**
* 2_1, 5_1, 7_1: |A| = |B| = 8. The 1–2 hidden lines have small κκ' (≈ 1e-4 in 2_1), so their deviations fall
  near or below the 1e-6 threshold (deg 5–7 instead of ≥ 50).

**All 28 completions fit the κ-family normal form** (residuals ≤ 2e-14). In every case A is a clique, every A×B pair
is non-orthogonal, and none is an A2 triple. By type:

| A, B type | completions |
|---|---|
| \|A\| = 8, B = A | 17 |
| \|A\| = 7, B = A | 6 |
| \|A\| = 6, \|B\| = 9 | 1 |
| \|A\| = 3, \|B\| = 19 | 1 |
| \|A\| = 1 (rigid side Y), \|B\| = 57 | 2 |

**Moduli tangent dimension = 2160 − rank J − 200** (GL9 × torus orbit, rank 200, all gaps ≥ 1e-6/1e-12). It **equals
the κ-family dimension |A| + |B| − 1 − [|A| = 1] at all 28 points** (13, 14, 15, 21, 56).

So every sampled point is a smooth point of its lift component, and no other family passes through it. This is
expSJ next-experiment 2, answered negatively: there is no tangent excess.

## 4. Non-clique supports (gap not in expSM §4) — CLOSED for independent A (`nonclique.py`, `nc2.py`, `ncctrl.py`)

The κ-form only needs A = S to have independent roots. A need not be a clique: lines of A∖B may be mutually
orthogonal, and only A∩B must be odd to all of A. expSM's orbit/Clean tables covered cliques only.

If A contains an orthogonal pair {s, t}, then B ⊆ C(A) ⊆ C2 = nonorth(s) ∩ nonorth(t) (24 lines).
* There is one W-orbit of orthogonal pairs.
* C2 has exactly 2 bad classes (486, 504). Both contain all of C2, so **non-Clean B forces B = C2 = C(A)**.
* No third line keeps C(A) = C2, so **A = {s, t}**.
* Every other non-clique independent A has Clean B, and Theorem 1 applies. Its proof uses only that the standard
  indices are pairwise odd, which holds.

For A = {s, t}, B = C2:
* The X side has only e9 (Def dim 7). g_s and g_t are not deficient, because s ⊥ t and α ≠ 0.
* If e9 is used, credit ≤ 7 + dim span(q_s, q_t) = 9. The other lines have γ = 0 and are non-deficient.
* If e9 is not used, T_Y ≤ 2 (e9') + 1 (the pencil line L_s or L_t) ≤ 3. Both pencil lines cannot occur together,
  because the ratio of signs is not constant.
* Exact control (`ncctrl.py`, 6 members, generic and special): PASS, mincost 72.

## 5. Status of C(P120) ≥ 72 after this run

**Proved (modulo the assumptions R1–R2 below).** Every square (m = 9) decoder of a rank-9 completion that is a lift on
a support S with independent roots has cost ≥ 72 in every GL9 gauge. This holds for every A, every B, and all special
parameters. Together with the A2-triple result (§2), it covers:
* every clique S (|S| = 1..8; the A2 triple is the only dependent clique, and no 9-clique exists);
* every non-clique independent S.

Sources:
* expSM Theorem 1 and (d);
* this run: §1 (|S| = 1), §2 (A2 non-concurrent), §4 (non-clique).

**Still open (exactly):**
* **R1. Classification.** Is every rank-9 masked completion on P120 a lift? The evidence is now stronger:
  * 28/28 saved completions are κ-lifts;
  * at all 28, the tangent dimension equals the family dimension, so each sampled lift component is locally the whole
    variety;
  * but components that contain no lift are not excluded.
* **R2. Rigidity.** Is every rank-8 completion on T = V∖S a root realisation? This is assumed by the Lift Lemma. It is
  proved for P112-type T only, and checked numerically for P120.
* **R3. Lifts on dependent non-clique S.** These are supports containing ≥ r+1 lines of a rank-r subsystem other than
  an A2 triple, e.g. 4 lines of A3, and every |S| ≥ 9. There, concurrency of the c_s can fail, as in §2. The normal
  form is not derived. None appeared in the census.
* **R4. Outside the rank-9 square family:**
  * rank-8 N with a 9-term decoder (case (b)) on P120. #2724 closes it on P112 and P104 only, and monotonicity goes
    the wrong way;
  * m ≥ 10 terms (#2338 §5.2(d); `freeh.py` is Boolean-SAT there, so only algebraic arguments help);
  * h ≥ 10.

**Next step (ranked).**
1. **R3 sweep.**
   * Enumerate W-orbits of dependent non-clique S with |S| ≤ 9.
   * For each, derive the c_s freedom exactly, as in `a2lib`: the kernel of the η-consistency system.
   * Where an extra parameter exists, rerun the `a2scan` pipeline: T-kernels → line lattice → `a2credit`-style
     enumeration.
   * Stop test: a credit-10 survivor, which is then built exactly in Fractions and checked with mincost.
2. **R1.** Run an LM census with forced non-lift ansätze, e.g. two rank-9 sides seeded off every lift component. Each
   new completion goes through `census3.py`.
   * Stop test: a completion that fails the κ-fit, or that has tangent > family dimension.
3. m = 10: lift the Theorem 1 label argument to rank-9 factor pairs with one extra rank-1 term.

## Files

| script | role | log / result |
|---|---|---|
| `s1lib.py` | \|S\| = 1 geometry: E, N, p_b, M_b^±, H0, h_b, C-parts, `dimD` | — |
| `geom1.py`, `geom2.py`, `facts.py` | Lemma A finite premises | `geom1.log` (its last line prints dim off by one: D(u0) = {0, H0, h_b, h_b*} is dim 2), `facts.log` |
| `s1relax.py [nob\|nodim\|norank]` | Case II CP-SAT | `s1relax.log`: 9; nob → 16 |
| `s1gs.py [ctrl\|noS\|sStar]` | Case III CP-SAT | `s1gs.log`: 9; noS and sStar → 10 |
| `hz.py` | H0-generic case | `hz.log`: T_Y ≤ 3 |
| `kctrl1.py seed n` | exact Lemma A control plus mincost | `kctrl1.log`: PASS 48, mincost 72 |
| `a2lib.py`, `a2scan.py`, `a2scan2.py` | A2 normal form, Lemma B | `a2scan*.log` |
| `a2credit.py [ctrl]` | A2 exact credit enumeration | `a2credit.log`: 9; ctrl → 11 |
| `a2ctrl.py seed n` | exact A2 members incl. special ones | `a2ctrl.log`: PASS 32, mincost 73 |
| `tangent.py`, `census3.py`, `hidden*.py` | lift fits plus tangent dimensions of all 28 | `census3.log`: all OK |
| `nonclique.py`, `nc2.py`, `ncctrl.py` | non-clique supports | `nonclique.log`, `nc2.log`, `ncctrl.log`: PASS |

All scripts are run as `../expO/run <script>`, single worker. CP-SAT runs take 40–70 s each, and `census3.py` takes ~4 min.
