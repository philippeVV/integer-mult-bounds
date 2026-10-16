# Recursion, tensor powers and cross-level overhead sharing in the local word: closed at ×1.000

This closure covers three routes, all ×1.000:
- **recursion:** a ZLL-style deep word;
- **tensor power:** an unrolled word of complete levels;
- **cross-level overhead sharing:** a k-level word paying the centre or helpers only once.

It merges three scouting runs: expQC, expQG (the primary source for the lemmas) and expQM (the loophole).

Status: scoped closure, author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`.
Multipliers are over κ = 0.0017. ZLL = Zhang–Li–Li, arXiv 2610.09367. It was read only as text; none of its code was run.

## Verdict

| Item | Status | Multiplier |
|---|---|---|
| **Tensor invariance**: Ψ_{S⊗S′} = Ψ_S·Ψ_{S′}, so a(S⊗S) = a(S) and a(S⊗S′) ∈ [min, max] | proved (new); checked exactly | **×1.000 per level** |
| **Lemma W**: s ≥ m·\|Fix ρ\| | the trivial-routing case of upstream `stronger-rank-screens.md` §2 (fixed strands form a unit integral routing of the fixed demands) | — |
| **Lemma W′ (new)**: s ≥ m·cyc(ρ) | proved over ℚ; tighter than W whenever ρ has cycles of length ≥ 2 | — |
| **Switch-cycle bound (new)**: s ≥ m·max_π cyc(ρ⁻¹π) over per-gate port permutations π | proved over ℚ. π = id gives W′; π = ρ gives the upstream routing theorem s ≥ Wm. **Vacuous on single-copy Eulerian schedules**: an Euler-circuit pairing gives cyc = 1. | — |
| ZLL's own circuit | integer, so `stronger-rank-screens` §1 gives s ≥ Wm | ×0 |
| ZLL seed, ignoring §1 | a ≤ ε₀/ln b₀, with log₁₀ ε₀ ≤ −1936 at the required h = 3 | ≈ ×10⁻¹⁹³⁰ |
| Compiled depth-j word S_j (compute–copy–uncompute) | δ_j ≤ k/W ≤ 1/(3 + 2(m₀/k₀)^{2^j}) by W′, which decays doubly exponentially | ×0 asymptotically |
| **Cross-level overhead sharing** (a k-level word pays the centre, or the centre and helpers, once) | Needs slot coincidence or frame coincidence. Slot coincidence is empty on #1094 and #589 (Kronecker overlap lemma plus an exact check). Frame coincidence is expNV's flat tower (n ≤ 4) and is batch-floored besides. | **×1.000** realised. Ideal ceilings: ×2.77 (centre once, k → ∞) and ×9.65 (centre + helpers once), both unreachable. Generous residual ≤ ×1.016. |

## Changed hypothesis tested

The local word is normally a constant-size seed compiled once. Here it is a **depth-j recursive network** (ZLL's S_j) whose
coding gap grows with j. The hope was that the gap grows faster than the supplier pays for it.

## Proof sketches

**Uniform-metric dictionary.** A uniform hop metric x_e = 1/b on a network with k demands, m edges and own-pair distance b
gives s/(Wm) = m/(kb) = ρ (ZLL's ceiling) and L_ent = ln b. Hence **a = ln(1/ρ)/ln b = 1 − ln(m/k)/ln b**.

**Tensor invariance (hypothesis H_T: complete supplier levels).** Let one level be a histogram (n_r, m, W) with
Ψ_S(a) = W⁻¹ Σ n_r (r/m)^{1−a}, and let a_S solve Ψ_S = 1.
- Unrolling two levels gives the histogram S⊗S′, with entries n_r n′_{r′} at rank rr′ and width mm′.
- Since (rr′/mm′)^{1−a} = (r/m)^{1−a}(r′/m′)^{1−a}, we get Ψ_{S⊗S′} = Ψ_S Ψ_{S′}. So Ψ^k = 1 ⇔ Ψ = 1.
- ln Ψ is convex and increasing, so a(S⊗S′) is a secant-weighted average of a_S and a_{S′}, which lies **between them**.
  Direct sums also lie between them.
- ZLL's ln-gap and ln-distance both double per tensor level, so a_j = ln(1/ρ₀)/ln b₀ for all j.

*Correction to the earlier scouting note.* That note stated a "mediant" formula Σ ln(1/(1−δ_k))/Σ ln m_k. It is exact only
when every level has a single rank. On random histograms it is off by up to 97% (`check_qg.py` Part 2). The correct statement
is "between a_S and a_{S′}", and the conclusion "≤ max a" is unaffected.

**Switch-cycle bound.** At every gate, choose a permutation of its ports, and let π be the composed role permutation. A port
permutation is a 0/1 scalar map, so it commutes with the gate's common frame. The upstream telescoping
(`stronger-rank-screens.md` §2) then gives s ≥ rank(D_out P_π − P_π D_in) = Σ_w rank(M_out(πw) − M_in(w)). Summed over a
cycle of σ = ρ⁻¹π of length ℓ, the endpoint identity turns the blocks into ℓ·I_m, which has rank m over ℚ. Each cycle
contributes ≥ m.

Hypotheses, each load-bearing:
- one role per segment;
- terminal segments charged;
- the cycle form needs char ∤ ℓ (the contract's frames are rational);
- an all-role permutation with one common frame per gate.

**Lemma W on ZLL compiles.** ρ is the identity except for k transpositions. So |Fix| = W − 2k and cyc = W − k, which gives
δ ≤ k/W = k/(3k + 2E). ZLL's tensor squares m/k per level.

## The loophole: pay a per-level overhead once (expQM), now closed

Invariance concerns *complete* levels. Now take a k-level word that pays an overhead X only once. Its levels are S ⊗ S_X^{k−1},
where S_X is S with X deleted. Then Ψ = Ψ_S Ψ_{S_X}^{k−1}, and the root a_k tends to a(S_X).

**Ledger: what could be paid once.**
- **In principle:** the centre/decoder copy charge (64 / 63 per call) and the helper climbs (65% of Λ on complex).
- **Never:** the data climbs (30% of Λ: each level climbs its own frames) and the unpaid coordinate (that is the gain).
- **Adds nothing:** the idle remainder, which is already ideal in the relaxation.

**Ideal prices (`price_qm.py`, pinned engines, both suppliers).** These are ceilings, unreachable as shown below.

| Word | k = 2 | k = 4 | k → ∞ | Reconciles with |
|---|---|---|---|---|
| centre once | ×1.89 | ×2.33 | **×2.77** | `levers.py:31` "centre path free" ×2.803 (complex side, to 4 digits); `levers.py:32` "shared by 2 calls" ×1.901 vs the 2-level root ×1.903 |
| centre + helpers once | ×2.93 | ×5.0 | **×9.65** | expM1 "every lever free" ×9.7 |

The joint ceiling is not ×23 = 2.8 × 8.2, because the data floor is paid at every level. Complex-only, it is ×1.148 against
the bit relaxation and ×1.018 against the realised bit endpoint.

**Kronecker overlap lemma (proved; `check_qm.py` P1 checks it exactly on 300 random instances over GF(10007)).**

  dim(A⊗F^{m′} + P⊗B) = dim A·m′ + dim P·dim B − dim(A∩P)·dim B,

because (A⊗F^{m′}) ∩ (P⊗B) = (A∩P)⊗B.

So one copy serving an outer centre (carrier slots A) and the centre B of an inner call hosted in the outer child P saves
exactly dim(A∩P)·dim B. That is zero unless the host child transforms the carrier's content on a carrier slot.

**No carrier/source overlap on the selected words (`check_qm.py` P2, P3, exact).**
- **#1094.** Every source register's climb into F_k precedes six gather adds, and its only same-content climb leaves F_k
  upward. Backward, the carrier is a fresh sum subtracted into sinks that already carry 12–38 earlier writes.
- **#589.** No retained total climbs between its last write and the snapshot.
- So A∩P = 0 for every host child outside the carrier subtree. Helper content is never copied, so helpers cannot be paid
  once either.

**Frame coincidence.** The other route is to re-embed the inner call into the outer frame chain. That makes it a flat word at
m = 73, which `../flat-tower` closes (n ≤ 4). Even granting a common chart, the inner inputs are independent of the outer
queries, so the p112 batch floor still charges 64 per demand.

**Generous residual.** Grant that every call hosted inside an outer centre child skips its own centre, at every depth. The
two-type branching moment then gives complex ×1.0158, bit ×1.0175, and **κ ×1.016**. No mechanism for it exists.

**Scope.** This uses the budget law's per-register child price. Interfaces that price a multi-register child as one transform
are outside it (budget law §3: ⊕ ≠ ⊗).

## Stop test that fired

The test was a(S_j)/a(S_0) ≤ 1 for all j, together with "the pre-filter bites no live candidate". Both hold:
- tensor invariance gives the first;
- the Lemma W/W′ pre-filter cannot bite on any helper-free single-copy Tanner schedule (Euler pairing: |Fix| = 0, cyc = 1);
- on the compiled candidates it gives δ ≤ 0.081 (BH) and 0.070 (PG(2,16)−M), where Lemma C (F* = W) binds anyway.

For the loophole, the test was: every shareable term needs slot or frame coincidence. It fired: slot coincidence is empty on
both selected words, and frame coincidence is the closed flat tower.

## Files and reproduction

Run each from this directory: `python3 -I -B <script>`. Everything is stdlib.

| Script | What it checks | Runtime | Output |
|---|---|---|---|
| `check_qg.py` | **Part 1:** W, W′ and the switch-cycle bound on 219 random permutation circuits scored by upstream `scripts/finite_bit_contract.py`, with a greedy adversarial frame minimiser and controls. **Part 2:** tensor invariance on the selected profile and 2,000 random pairs, and the mediant counterexample. **Part 3:** the pre-filter on `inputs/`. | 10 s | `check_qg.out`, `PASS_EXPQG` |
| `zll_scale.py` | the HWZ/ZLL ratio recursion, ε₀ and Lemma W ceilings, the S_j table, and the uniform-metric dictionary on BH and PG(2,16)−M | 0.1 s | `zll_scale.out`, `PASS_EXPQC` |
| `check_qm.py` | P1: Kronecker overlap lemma. P2: carrier/source overlap on the real #1094 programs (`research/community-frontier1094-audit/selected-replay.tar.gz`). P3: the same on the real #589 word (`research/e8-bit-p111-tuned/tune/`). | 2.5 s | `check_qm.out`, `PASS_EXPQM_CHECK` |
| `price_qm.py` | ideal amortisation (centre once; centre + helpers once) and the generous residual, via `engines.py` (the pinned `ceiling.py` and `levers.py`) | 108 s | `price_qm.out`, `PASS_EXPQM_PRICE` |

`check_qg.py` imports `scripts/finite_bit_contract.py` and `scripts/audit_joint_frames.py` from the repository root, found
relative to this file. Override it with `--repo PATH` or `IMB_REPO`.

`inputs/` (≈ 150 KB) holds the candidate schedules, read-only:
- `schedule784.json` (expNI);
- the five expNC causal/PSA schedules, converted from pickle to JSON with keys `pos`, `R` and `edges`;
- `bh_compile.json` (expNN).

The outputs match the scouting logs except for the candidate names, which lost their `.pkl` suffix.

## Overlap

`find_overlap.py` and `gh search` returned 0 hits for tensor invariance, recursion-depth invariance, fixed-point full width,
and own-strand telescoping. The nearest items are:
- upstream `stronger-rank-screens.md` §2, of which Lemma W is the trivial-routing case;
- #1128 (full-width endings for the #589 helpers);
- #2429 (T74 restoration climbs).

If this is posted, it belongs as a comment on the budget-law thread (#2278), not as a standalone PR.

Prepared with Claude Code (Claude Opus 5.5) assistance. Not independently reviewed.
