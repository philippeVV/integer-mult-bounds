# Orthogonal-frame plane upgrade (partial u = 2): fails in the natural decoder classes

**Scoped obstruction; no κ claim.** This tests the experiment proposed in discussion #2278. Upgrading a
set S of ports from lines to planes ⟨q_a, r⟩ raises κ iff the extra decoder charge is
ΔC < (8/9)|S|. Here S is an orthogonal 8-frame with a shared r, so the bar is ΔC ≤ 7.

Pricing follows #446 Theorem C / #1516: a term d tᵀ costs the F2 pairing rank of the plane spans
of its supports. Within a frame with shared r, two planes pair as I₂, so any term touching S×S costs at least 2.

`check_frames.py` is standard-library only and prints `PASS_FRAME_UPGRADE_CHECK` (output in `check_frames.out`).
It checks every claim below on all 1260 frames × 16 valid r = 20160 cases from the selected labels.

## Results
1. **Structure of every case.**
   - Size of T and lost entries per side: |T| = 48 with 192 lost, |T| = 52 with 208, or |T| = 56 with 224.
   - Cross-block ranks: rank_Q G[S,T] = rank_Q G[T,S] = 8, and the F2 rank is 4.
   - Every b ∈ T loses exactly 4 frame ports.
2. **Minimum rank is exactly 16.** T contains a pairwise-even 8-set S″, so U = S ∪ S″ has
   N[U,U] = I₁₆ and α(Γ′) = 16 in all cases. The decoder I_S ⊕ (G/2)[S^c] has rank 16.
3. **Floors.**
   - **Any decoder:** at least 16 terms of positive cost.
   - **Core class**, where some 8 terms are worth at least 64 in total (e.g. any re-choice of the
     incumbent hyperplane terms): **ΔC ≥ 8 > 64/9.**
   - **Additive class**, a rank-8 gauge of the incumbent plus corrections: the forced blocks
     X[S,T] = X[T,S] = −G/2 give rank ≥ 16, so **ΔC ≥ 16**. This is attained.
4. **Best exact construction:** I_S ⊕ (G/2)[S^c], cost 80 (ΔC = 16), exact on 47 random cases.
   Annealing over class decoders never beat ΔC = 2|S|.
5. **Break-even** (`be.py`, pinned `ceiling.py`): for |S| = 8, ΔC = 7 wins by +0.18% and ΔC = 8 loses.

## Still open
- A decoder with at least 16 terms averaging ρ ≤ 4.4, whose eight costliest terms total at most 63.
  This is outside both classes and is the same low-ρ regime as the general problem.
- Other choices of S that the rank argument does not exclude in the core class (numerical
  evidence only):
  - Hamming classes with shared r: α = 15 in 56 of 224 cases. A rank-15 decoder exists
    numerically, leaving a hairline +0.18% case.
  - Two frames, each with its own r.
  - Per-port r_a.
  All of these fail in the additive class, and the best construction is ΔC = 2|S|.

Prepared by Philippe Van Velzen with substantial Claude Code (Claude Opus 5.5) assistance;
no independent human review is claimed. Apache-2.0.
