# κ and undirected coding gaps: oblivious multipliers, the AFKL transfer, and the BH/BGS exponent

FINDING / THEORY NOTE. No κ claim. Author-checked, not independently reviewed.

1. **Theorem A (`AFKL-TRANSFER.md`).** The network-coding conjecture enters the multiplication argument of Afshani–Freksen–Kamma–Larsen (2019) exactly once, and linearly. Replacing it with "flow ≥ coding/γ" gives a bound for every fan-in-c circuit computing the n-bit shift: s ≥ n·log_{2c} n / (4γ(s)).
2. **Corollary B.** An oblivious multitape TM running in O(n lg^{1−κ} n) yields an undirected coding/flow gap family with exponent η ≥ κ. Equivalently, κ_obl ≤ η_* ≤ 1. So κ → 1 *requires* coding gaps that saturate the flow–cut bound. This is a necessary condition only.
3. **Obliviousness (`HVDH-OBLIVIOUS-AUDIT.md`).**
   - Upstream's multiplier and Harvey–van der Hoeven 2021 §§2–5 are oblivious after fixed-width padding at ×O(1) cost. The data-dependent spots are:
     - variable-length integers;
     - the Lemma 2.2 scaling branch (fixed by clamp-and-select);
     - sign handling.
   - This closes the gap #634 left open ("HvdH internals not established").
   - Consequence: if upstream's Theorem 1 holds at κ, then γ(N) = Ω((lg N)^κ) at every N.
4. **First explicit exponent for the Braverman–He seed amplified by BGS.**
   - c₁ = 67.94, giving Theorem 1.1's c₂ = 0.00298.
   - The all-N liminf is 0.00198, or 0.00260 with α optimised.
   - Today's κ ≈ 0.0017 is 0.66–0.87× of the all-N guarantee. So **×1.16–1.53 in κ would exceed the best published all-N coding-gap guarantee from that seed.**

Reproduce (standard library, each under 1 s):
- `python3 exponents.py`
- `python3 bh_c1.py`
- `python3 bgs_dense.py`

Credits: #634 (ShikeChen01) for Lemma B′, the HvdH Lemmas 2.1–2.5 audit and the "κ ≤ γ" question. Papers: AFKL (arXiv 1902.10935), Braverman–Garg–Schvartzman (arXiv 1608.06545), Braverman–He (arXiv 2610.10108), Zhang–Li–Li (arXiv 2610.09367), Harvey–van der Hoeven (Annals 2021).

Prepared with Claude Code (Claude Opus 5.5) assistance, using sub-agents; the main session re-ran the numeric scripts. Not independently reviewed.
