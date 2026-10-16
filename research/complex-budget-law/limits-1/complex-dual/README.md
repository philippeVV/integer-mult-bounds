# Complex dual (expXG, expXK): a gauge-invariant ϑ-type bound for Problem X(ℂ), the proved ceiling map by girth, and Conjecture Θ

Status: scoped ceiling work, author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`. Multipliers are
over κ = 0.0017. Extends `../../closures-3/complex-coding-seeds/` (expRF T1–T3, which covered flat-leverage codes only).

## Verdict

| Item | Value |
|---|---|
| Question | is s = dim W/E ≤ 1/g (equivalently a ≤ 0, with a = 1 − ln((1−s)/s)/ln(g−1)) for **every** MDS half-rate complex Tanner code, not just flat-leverage ones? |
| Gauge (Lemma D) | dim W is invariant under positive diagonal rescaling D of the edges; leverages and moments are not. Along D_t = t^{1_{X^c}} (X an n/2-factor) **E/2 eigenvalues of K = 2I − Π_P − Π_L go to 0 for every code**. Any valid dual must take an inf over gauges: Θ̂* = inf_D Θ̂_J(D·S) |
| Girth 4 (Thm B) | **closed for all codes**: Barthe/Forster scaling flattens one side, then T1 gives s ≤ ¼ = 1/g. Tight on K₄,₄ |
| Any field (Thm C) | the gauge boundary is the bootstrap bound dim W ≤ Σ_v (deg v/2 − #earlier nbrs)^+, valid in every characteristic (`../bootstrap/`) |
| Flat codes (Thm F) | the kernel atom ν({0}) ≥ ½ refines the Markov bound to **s ≤ 1/g exactly at g = 4, 6, …, 16**; expRF's girth-6/10 slivers (2/11, 0.1034) close |
| Girth 8 (expXK) | **two-sided flat ⇒ Θ̂₃ = 1/8 for any codes** (Thm 2, via the exact m₃ identity, Lemma 1); the complementary gauge exists (Thm 3) but cannot prove the general case (Prop 4); ETF codes attain 1/8 (Prop 5) and are local maxima |
| Conjecture Θ | for g ≡ 0 mod 4, every MDS code has a gauge with Θ̂_{g/2−1} ≤ 1/g. **Open at girth 8, and tight** (attained by SIC/icosahedral ETF codes; no adversary exceeded 1/8 on W(3), W(4)−M, W(5)) |
| Girth 6, 10 | **open for non-flat codes**, and the local even-J dual is provably insufficient at g = 6 (codes keep Θ̂ > 1/6 on every bounded gauge box) |
| Multiplier | **×0** realised (ceiling work) |

## Proved ceiling map for the complex supplier (a ≤ 0 ⇔ s ≤ 1/g)

| girth | flat / transitive codes | codes with a two-sided flat gauge | all MDS complex codes | every field (bootstrap, graph-specific) |
|---|---|---|---|---|
| 4 | closed (T1) | closed | **closed** (Thm B) | — |
| 6 | **closed, s ≤ 1/6** (Thm F) | **closed** (Thm F: m₁ = ¼, m₂ = 3/16 for any two-sided flat codes) | **open**; the moment dual alone cannot close it | closed on PG(2,3), biaffine(4), PG(2,4)−M, PG(2,5) (Boot/E ≤ 1/6) |
| 8 | closed (Thm F) | **closed** (Thm 2: Θ̂₃ = 1/8) | **open = Conjecture Θ** (tight; evidence on W(5), the first graph with Boot/E > 1/8) | closed on W(3) (= 1/8), W(4)−M |
| 10 | closed (Thm F, 1/10) | not extended (needs J = 4) | open | none tested |
| 12, 16 | closed (Thm F) | not extended (J = 5, 7) | open | — |

For κ ≥ 0.1 the target is a ≥ 1/9, not a > 0, so every "closed" cell above also closes κ ≥ 0.1; the open cells are where
a complex Problem X solution could still live. Constructive hint from the dual: at girth 6 the maximising profile is
non-flat with positively correlated leverages (m₁ ≈ 0.244 < ¼, m₂ ≈ 0.187 > 3/16), on graphs with Boot/E > 1/g.

## Changed hypothesis

expRF bounded only codes with constant leverage. The question was whether a Lovász-ϑ-type relaxation (moments of
Π_PΠ_L, Markov/Christoffel certificates) bounds all codes. It does, once the dual is made gauge invariant, at girth 4 and
for two-sided-flat codes; elsewhere it stops exactly at the gauge boundary, which is the field-free bootstrap bound.

## Definitions and proof sketches (expXG)

**Θ̂.** m_k = tr((Π_PΠ_L)^k)/E for k ≤ J = ⌊(g−2)/2⌋ (closed alternating walks are tree-like), m′_k = 2m_k. Then
Θ̂_J = ½·max{μ′({1}) : μ′ on [0,1] with moments m′₀..m′_J} = ½/(vᵀH⁻¹v), with v = 1 and H the Hankel matrix (Lukács
certificates p = q², x·q², q(1) = 1). **Validity:** W is the eigenvalue-1 space of Π_PΠ_LΠ_P, whose rank ≤ E/2 forces
ν({0}) ≥ ½, so s = ½ν′({1}) ≤ Θ̂_J for every D. ∎

**Theorem B.** MDS ⇒ ½·1 ∈ relint of the matroid base polytope ⇒ a positive scaling with all point leverages ½ (Barthe
1998, Forster). Point stars partition E, so all points flatten at once; T1 gives s ≤ Σ a_e b_e/E = ½Σ b_e/E = ¼. ∎

**Theorem C (any field).** At v with j back-edges T, extend T to an information set I of A_v; W|_{E(v)} embeds in F^I, so
the new span at v is ≤ (k − j)^+. The same bound is the limit of Θ̂* along degenerate gauges t^c. ∎

**Lemma D.** D_t sends every A_v to F^{X∩E(v)}; K(D_t) → diag(0 on X, 2 off X); |X| = E/2 eigenvalues → 0 while dim W is
fixed. ∎ (So the expXE r ≥ 9 "solutions" on PG(2,3) were gauge artefacts; Thm C gives dim W ≤ 8 there, and expXE's r = 8
code is a genuine MDS optimum, `xg/check_xg.py` C3.)

**Theorem E.** (A⊗A′) ∩ (B⊗B′) = (A∩B)⊗(A′∩B′); product certificates and gauges tensor, so Θ̂(S⊗S′) ≤ Θ̂(S)Θ̂(S′). ∎

**Theorem F.** For the free law ν′ is arcsine, m′_k = C(2k,k)/4^k; its degree-J Christoffel mass at 1 is 1/(J+1), so
Θ̂ = 1/(2J+2) = 1/g. With both sides flat and girth ≥ 6, the tree expansion gives m₁ = ¼ and
E·m₂ = Σ_p b_pᵀ|Π_p|^{∘2}b_p + Σ_l a_lᵀ|Π_l|^{∘2}a_l − Σ_e a_e²b_e² = E/8 + E/8 − E/16 = 3E/16 for any codes. ∎

## Girth 8 (expXK)

**Lemma 1 (girth ≥ 8).** With D_v the opposite leverages on E(v), K_v = Π_vD_vΠ_v,
u_e = Σ_{e′∈E(p)∖e}|A_{ee′}|²b_{e′}, w_e = Σ_{e′∈E(l)∖e}|B_{ee′}|²a_{e′}:
E m₁ = Σ_e a_eb_e; E m₂ = Σ_v tr(Π_vD_v)² − Σ_e a²b²; **E m₃ = Σ_v tr(Π_vD_v)³ − Σ_e a³b³ + 3Σ_e u_ew_e.**
*Proof.* Closed alternating walks of length ≤ 6 are tree-like; the only side-mixing walk is centred at an edge, has 3
rotations and weight |A_{e₀e₁}|²|B_{e₁e₂}|²a_{e₂}b_{e₀}. ∎ Hence the local identity, with F = 4m₃ − 4m₂ + m₁ ≥ Θ̂₃ ≥ s:

  E·F/4 = Σ_v tr K_v(K_v − ½)² − Σ_e a_eb_e(a_eb_e − ½)² + 3Σ_e u_ew_e.  (★)

**Theorem 2 (two-sided flat ⇒ 1/8).** a ≡ b ≡ ½ gives K_v = ½Π_v (first sum 0), −1/64 per edge, u = w = 1/8 (+3/64 per
edge): F = 1/8 and m = (¼, 3/16, 5/32), the arcsine moments, so Θ̂₃ = 1/8 for any codes. ∎

**Theorem 3 (complementary gauge).** Φ(x) = Σ_v log det(Q_vᴴe^{2x}Q_v) − 2Σ_e x_e is convex (Cauchy–Binet),
∂Φ/∂x_e = 2(a_e + b_e) − 2, and attains its minimum because ½·1 ∈ relint B(U_{k,2k}). So a gauge D_c with a_e + b_e = 1
exists (E equations, E unknowns). At D_c, (★) is a sum of single-vertex terms plus edge couplings. ∎

**Prop 4 (D_c is not enough).** Second order at D_c: E(F − 1/8)/4 = δᵀMδ + O(|δ|³) with
M = −(3/16)I + ½(O_P + O_L) − (3/2)(O_PO_L + O_LO_P); near reducible (pairing) codes M|C has λ_max = +0.032 and
Θ̂₃(D_c) = 0.125001 > 1/8, while Θ̂* ≈ 0.08 there (reached by the degenerate gauge). Any proof must mix D_c with gauge
degeneration (bootstrap structure).

**Prop 5 (gauge Hessian at flat codes).** Θ̂₃(D_y) − 1/8 = (4/E)yᵀΓy + O(y³), Γ explicit (expXK §2). Fourier and random flat
codes: Γ indefinite, so nearby codes have Θ̂* < 1/8. ETF codes (qubit SIC on W(3), icosahedral ETF(3,6) on W(5)): Γ ≻ 0 on
1^⊥ and Θ̂* = 1/8 exactly. At SIC codes the code/gauge Schur complement has 559 zero modes (400 flat-manifold tangent + 159
gauge orbits) and 81 strictly negative normal directions (`xk/sic_schur.py`): **Conjecture Θ holds to second order near its
maximisers.**

**Adversary (numerical).** max over codes of inf over gauges, J = 3: W(4)−M ≤ 0.1031; W(5) (degree 6, E = 936,
Boot/E = 0.137 > 1/8) 0.12499917 from random, 0.124986 from icosahedral starts; W(3) from SIC 0.12499993. No code with
Θ̂* > 1/8. At girth 6: PG(2,3) codes keep Θ̂ = 0.1723 > 1/6 on the gauge box |log d| ≤ 1 (0.1671 at 1.5); the unconstrained
inf runs to the bootstrap value 8/52 = 0.1538.

## Failed proof routes (recorded so they are not repeated)

Per-vertex LP splitting of (★) (n = 6 LP value 2.1e-4 > 0; Prop 4 explains why); AM–GM 3uw ≤ (3/2)(u² + w²) (fails at
vertex level); flattening both sides at once (overdetermined, V(n−1) > E); "first-order gauge optimality forces flat
leverage" (false: interior optima keep leverages in [0.27, 0.76] on W(3)).

## Next step and stop test

1. Prove Γ ⪰ 0 ⇒ normal Schur complement ≺ 0 at every flat degree-4 code (local Θ). **Stop/refute:** a flat degree-4 code
   with Γ ⪰ 0 and a positive normal Schur eigenvalue refutes Θ locally; certify its dim W.
2. Test the mixed gauge inf over {D_c·t^c} (vertex orders c) on the adversary end points and the Prop 4 codes; if always
   ≤ 1/8, prove the inequality from (★) plus Theorem C.
3. Girth 6/10: no dual-only route exists; any proof needs rank/bootstrap structure. Constructions (if attempted) must be
   non-flat, positively correlated, on graphs with Boot/E > 1/g (PG(2,q), q ≥ 17; W(q), q ≥ 11). `../seeds-negative/` (expXM)
   records that Singer-equivariant codes on PG(2,17) fail by ×8.4.

## Files and reproduction

Runtimes under a 4G cap, one worker.

| Script (run from its subdirectory) | Checks | Runtime | Output |
|---|---|---|---|
| `xg/`: `uv run -q --no-project --python 3.12 --with numpy --with scipy python -I -B check_xg.py` | 29 checks: C1 free-law values 1/g (g = 4…16) and Fourier codes; C2 Theorem B on PG3/BA4/W3 (flattening, Σa_eb_e = E/4) and K₄,₄ tightness; C3 bootstrap closures (PG3, BA4, PG4M, W3, W4M, PG5) and the expXE r = 8 MDS optimum (dim W = 8, gap 0.293); C4 Lemma D; C5 tensor product | 3 s | `check_xg.out`, `PASS_EXPXG` |
| `xk/`: `uv run … --with numpy --with scipy python -I -B check_xk.py` | 14 checks K1–K10: sparse = dense Θ̂₃; Lemma 1 against exact traces (W3, W5); Thm 2; Thm 3 and (★); W5 adversary end point < 1/8; analytic Γ; ETF codes Θ̂* = 1/8; Prop 4 | 18 s | `check_xk.out`, `PASS_EXPXK` |
| `xk/`: `uv run … --with numpy --with scipy python -I -B sic_schur.py` | Schur complement at SIC codes on W(3): 559 zero modes = 400 + 159, 81 negative normal directions | 12 s | `sic_schur.out`, `PASS_EXPXK_SCHUR` |

Included but not replayed (scouting): `xg/xg_core.py` (Θ̂, gradients, Barthe scaling; library), `xg/boot.py`,
`xg/artifact*.py` (gauge-artefact check of the saved expXE r = 9 point, `sol_PG3_r9_s0.npy`), `xg/t3.py`; adversary logs
`xg/adv*.log`, `xk/adv_*.log`; `xk/outer_xk.py` (adversary), `comp.py`, `local.py`, `gauge_hess*.py`, `maxmin_hess.py`,
`sic.py`, `mk_init.py`, `hess*.py`, `probe.py`, `coupled*.py`, `evalbest.py`, `etf_w5.py`, `boot_w5.py`, the failed LP route
`split*.py`, `vertex.py`, `verify_cert.py`, `flat_hess.py` with their logs and `.npy` certificates. Graphs: `xg/graphs_small.json`
(PG3, W3, BA4, PG4M, W4M, PG5), `xk/graphs_xk.json` (W3, W5, W7; copy of expXE `graphs_gq.json`).

Sources: expXG, expXK (expXE's r = 8 solution is copied as `xg/sol_PG3_r8_s0.npy`). Prepared with Claude Code (Claude Opus
5.5) assistance. Theorems B, C, E, F, Lemma 1, Theorems 2, 3 are elementary but not independently reviewed.
