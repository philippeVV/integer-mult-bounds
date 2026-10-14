# Plane-port decoder screens (follow-up to Discussion #1587)

**Scoped partial results. No κ claim; the question stays open.**

In #1587, ctoth's team showed that a plane-port local contract leaves two unpaid
coordinates per live family. Port t gets a plane `E_t = ⟨q_t, r_t⟩ ⊂ F2^9`; its
source runs `E_t → FULL` and its target runs `0 → E_t⊥`. The coverage identity then
becomes `mass = 71W + C·N`, and a decoder beats the incumbent iff its copy charge
`C < 1472/9 ≈ 163.6` per call. Under `ceiling.py` with `unpaid=2`, the largest
integer C that still wins is 160–162. This package screens both directions.

## Model

`lower/MODEL.md` sets out the model.
- The decoder needs `N = Σ d_k t_kᵀ` with `N_tt = 1` and `N_st = 0` whenever
  `E_s` is not orthogonal to `E_t`.
- Plane-orthogonal pairs are free and form a subset of the odd root pairs.
- A term costs the F2 pairing rank between the span of its row planes and the
  span of its column planes (#446 Theorem C, #1516).
- The 112 ports are identified exactly with the P112 root lines: all 14,400
  pair parities match.

## Results

**Proven or exactly checked:**
1. Each port keeps at most 28 of its 52 odd neighbours as plane-orthogonal,
   over all 128 nondegenerate planes through `q_t`, for all 112 ports
   (`construct/perport.py`, `lower/verify.py`).
2. **Floor.** `C ≥ α(Γ_E) + 1`, where Γ_E is the plane-orthogonality graph.
   The best plane choices found have α = 31, so this floor is at most about 32
   and **does not kill the route**.
3. **The tensor route is impossible.** A tensor decoder `d_t = c_t ⊗ a_t` with 2
   classes would give C ≤ 144. It needs a 2-colouring of the P112 odd graph with
   no monochromatic 5-clique, because five pairwise-orthogonal nondegenerate
   planes would need 10 > 9 dimensions. The odd graph has 161,280 5-cliques,
   and SAT reports that no such colouring exists.
   `construct/k5_independent.py` rebuilds E8 and the cliques from scratch and
   re-solves; it needs `python-sat`.
4. **Frame plus one star correction per port costs at least 176.** Every port
   has a broken row, and König's theorem forces at least 112 corrections on
   top of the 64 for the frame terms. Built exactly, it costs 288.
5. Keeping the incumbent E8 decoder only on Γ_E gives an N of full rank 112.
   The literal construction costs at least 184.

**Constructions:** every decoder built and checked exactly on the full instance
costs at least 224, the naive cost of one private rank-2 term per port.

**Numerical only (floating-point, uncertified). Corrected; see the third follow-up below.**
- The earlier ALS estimate of 56–90 for the minimum rank of N was an overestimate.
- Levenberg–Marquardt (`construct2/lm.py`) finds a symmetric rank-36 completion on the rW instance, with residual about 4e-22, and fails at rank 31.
- These completions are dense (ρ ≈ 9 per term), so they do not give a cheap decoder.

## Gain analysis (`gain/`, follow-up)

**Proven, with exact checks:**
6. **Theorem A (`integral_bound.py`).** Suppose every term `d_k t_kᵀ` is 2-integral, meaning no even denominators. Then `C ≥ 224`.
   - Reduce mod 2: `N mod 2 = I` on required entries, and the Gram `B` vanishes on free pairs.
   - So `Σρ_k ≥ rank_F2(B ∘ (I ⊗ J)) = rank blockdiag(Gram E_a) = 224`.
   - Subadditivity lemma: `Σρ_k ≥ rank_F2(B ∘ (Σ α_k β_kᵀ ⊗ J))` for any `α_k ⊆ supp d_k` and `β_k ⊆ supp t_k`. It is also checked on random term families.
7. **Gain inequality (`gain.py`).** `C ≥ 224 − gain(K₋)`, where `K₋` is the set of terms with even denominators. Their gain is (the F2 rank they remove) minus (their cost).
   - **A sub-160 decoder needs gain > 64.** Rebuilt exactly, the line incumbent has gain 48 (8 terms, cost 64, mod-2 pattern of rank 112).
8. **Top-layer lemma.** Let m be the largest 2-adic exponent.
   - The leading layer `S_m = Σ_{e_k=m} d̄_k t̄_kᵀ` vanishes mod 2 on every required entry.
   - Hence any split of that layer into two parts gives equal masked ranks, and the layer costs `≥ 2·max rank`.
   - Every diagonal covered at the top layer is covered by at least two top-layer terms.
9. **Why lines reach 48 (`ham.py`).**
   - Each port lies in exactly 4 frame terms. Those support sets are the 14 weight-4 words of the extended Hamming [8,4,4] code, with 8 ports each.
   - All overlaps are even, so the top layer cancels.
   - Each term's support is a whole hyperplane, which costs 8.
10. **Frame skeleton plus half-integer corrections costs ≥ 116 (`corrbound.py`).**
    - Choose each `r_a` inside the 4 frame hyperplanes that contain `q_a`. Then the frame still costs 64, but 2016 odd entries are lost.
    - Correcting them costs ≥ 52 by Theorem A applied to `2·N_corr`. The best exact build is 312 (`rectcover2.py`).
11. **Rank-≤2 Hadamard decoders `N = (G/2)∘M` are impossible for every plane choice (`k5.py`).**
    - The lost-pair graph would have to be bipartite, so the odd graph would need a 2-colouring with no monochromatic K5.
    - This is the same UNSAT as result 3, re-checked independently.
12. **Hamming classes need no decoder of their own (`cls14.py`).** With `r` inside the frame hyperplanes, every class pair has pairing rank 2, so each class costs exactly 16 (its naive cost).
13. **Floor `C ≥ |S| + |Q|` (`floor2.py`).**
    - S is a set of pairwise non-orthogonal ports, and Q ⊆ S is a set with pairwise pairing rank 2.
    - The floor is 103–108 on edge-maximising planes but only about 40 on α-minimising planes, so it does not close the gap.

**Searched, nothing found:**
- 2-term polarisation gadgets (`polar.py`) have gain ≤ 0.
- Simulated annealing over class-Hadamard decoders always collapses to 224.
- Frame-stabiliser-equivariant plane choices: none exist. The stabiliser has order 344,064, and we also tried all vector stabilisers and 109 vector-pair stabilisers. This was a sympy screen and is not packaged.
- Numerically, the masked rank on the frame-hyperplane instance is about 38–40, against α(Γ_E) = 31 (CP-SAT optimal).

**Gain status:** not settled. The route looks dead, but this is not proven:
- No construction beats gain 0.
- A general proof of `gain ≤ 64` would need to control the carry layers below the top 2-adic layer.

The scripts in `gain/` read the ports through `lower/ports.py`. That file loads the selected canonical word from `research/community-frontier1094-audit/selected-replay.tar.gz` and checks its hash. `k5.py` needs `python-sat`; `corr_rank.py` needs numpy.

## 2-adic Smith lemma (`smith/`, second follow-up)

**Proven for every plane assignment.**
- Let r_j(X) be the number of Smith invariants of X over Z_(2) with valuation ≤ j. Write e_k = v₂(d_k) + v₂(t_k), and let B̃ be the 0/1 Gram of the plane bases.
- Then `224 ≤ Σ_k r_{max(0,−e_k)}(B̃[R_k, A_k])`.
- Proof sketch: `B̃∘(2^m N⊗J) = 2^m I₂₂₄` holds exactly, and r_j is subadditive.
- At depth 0 this recovers Theorem A. The gain is bounded by Σ(h_k − ρ_k), where h_k is a per-term cap.
- Cap at depth d: `r_d ≤ Σ_{s=1}^{d+1} C(ρ,s)`. Exact caps over subspace pairs are 37/93/163 for ρ = 9 and 29/64/99 for ρ = 8, at d = 1/2/3.

| decoder class (any plane assignment) | C ≥ |
|---|---|
| all terms half-integral | 55 |
| half-integral, ρ ≤ 3 | 112 |
| every non-2-integral term has ρ ≤ 2 / ≤ 3 / ≤ 5 | 150 / 96 / 37 |

**Other results:**
- **Universal α(Γ_E) ≥ 14**, so C ≥ 15 for every assignment (`lagr.py`).
  - Average over the 2295 Lagrangians of j^⊥; every nonzero vector lies in exactly 135 of them.
  - The Hoffman bound on all 16,320 planes caps this averaging method at 15.
- **Per-instance floor** `C ≥ α(Γ_E) + α(Γ_{≤1})` (`floor3.py`). For example, rW gives 31 + 15 = 46.

**Still open.** Only decoders with a term of denominator ≥ 4 and ρ ≥ 7, or denominator ≥ 8 and ρ = 6, escape every universal bound above 32. That is exactly the line incumbent's mechanism: quarter terms of ρ = 8, with r₂ = 22 and Smith sum 176 ≥ 112 (`line_check.py`).

**Run** (Python 3.12 with numpy, scipy, ortools and python-sat):
- `smith/smith_check.py`
- `smith/corollary.py`
- `smith/lagr.py`
- `smith/line_check.py`
- `smith/floor3.py gain/rW_1.json`

## Third follow-up: large-term obstruction (`construct2/`)

- **Large terms need few independent ports.** Every term has rank 1, so K ≥ rank N ≥ α(Γ_E).
  - If every term has ρ ≥ 7, then C ≥ 7·α(Γ_E). A decoder below 162 built from large terms therefore needs a plane choice with **α(Γ_E) ≤ 23**.
  - For symmetric PSD decoders (N = Σ c_k c_kᵀ/2^j), Lovász gives K ≥ ϑ(Γ_E). For the line case this is tight: ϑ = 8 = rank G/2.
- **The best plane choice found has α = 28** (`construct2/amin_s1.json`, certified optimal by CP-SAT) and ϑ ≈ 32.5.
  - ϑ is computed by `theta.py` and is numerical, with no dual certificate.
  - On that instance, large-term decoders cost ≥ 196, and PSD decoders ≥ 231.
  - rW has α = 31 and ϑ ≈ 32.4.
- **Structured searches found nothing below 224:** the two-copy root decoder (G_ρ + G_σ)/4, hyperplane-supported quarter terms, and joint ALS.
- **Update: α = 25 is reachable** (`construct2/amin25.json`). CP-SAT proves this optimal; the graph has 756 edges and exactly 6 maximum independent sets, and ϑ ≈ 29.6. The search found it with tabu moves over single ports and port pairs and then stalled at 25. Large-term decoders on this instance still cost at least 7·25 = 175.
- On the lower-bound side, Lagrangian averaging is exactly 14 even with integrality. A generalised pairwise-dot-1 family was tried with CEGAR, with no result beyond 14.
- **The deciding question:** does any plane choice have α(Γ_E) ≤ 23? Averaging over Lagrangians only proves α ≥ 14, and real instances reach only 21–23 per Lagrangian clique.

## Open

**Decoder target: resolved in favour of the masked target (`target/line_decoder.py`).**
- An exact ℚ replay of the selected canonical word splits `y += Jx` at the cut into a drop-carried part N, through the 8 centres, and a lossless part L. It gives:
  - `J = I`;
  - `N = G/2`: diagonal 1, **0 on all 6,608 even pairs**, ±1/2 on odd pairs;
  - L supported on label-orthogonal pairs only, as Theorem A requires.
- So the incumbent decoder realises the masked identity, not the label Gram. `B₀ = (J + 2N) mod 2` is only the PR717 label-form identity.
- PR #1825's F2 factorization of `B₀` (charge 126) therefore does not meet the scalar endpoint. Open PR #1889 audits the same failure independently.
- All bounds in this package use this target.

The gap between 32 and 224 is open. Still unexplored:
- correction terms whose rectangles cover several rows;
- `r_t` taken from exterior coordinates;
- multi-total copy trees;
- other punctures.

A kill would need a bound on cancellation between cheap low-span terms.

Credit: ctoth's #1587 team posed the question. Background comes from #1417
(aliromman), #1516, #1348, #446/#595 and #968.

Prepared by Philippe Van Velzen with substantial Claude (Anthropic) assistance.
Two parallel sub-agents did the work; the main session re-checked results 1 and
3 independently. Apache-2.0.
