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

**Numerical only (floating-point, uncertified):** the minimum rank of N looks
like 56–90. A below-break-even decoder would therefore need many terms of
average cost about 2–3.

## Open

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
