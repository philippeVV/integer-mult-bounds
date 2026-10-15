# expSB: a second unpaid coordinate per family without generic plane ports

Scout for the "two unpaid coordinates per family" lever (levers.py ×3.77). Base: main `0fdd7d82`, with Discussion #1587 re-read in full on 2026-10-15 (`d1587-new.txt`).
Status: local only, no κ claim. `verify.py` checks every exact statement below and prints `PASS_EXPSB`. `lm_results.txt` holds the floating-point screens, which are evidence only.

## Changed hypothesis (one paragraph)
Today each live family has exactly one unpaid coordinate, the port label q. The centre (8 rank-8 copies, C = 64) serves the 112 unpaid port coordinates by realising the masked identity N = G/2. I tested four ways to give families a second unpaid coordinate without a generic plane-port label:
- (a) pair families so that they share it;
- (b) change the label alphabet (F4 or F8 lines instead of F2 lines);
- (c) take the coordinate from the recursion's output side;
- (d) upgrade only some families.

The existing plane-port obstructions (Theorem A, the Smith lemma, rank ≥ α+2, the 128 barrier) are proved for one plane per port. They do not directly cover shared planes, 3-dimensional F8 labels or mixed label dimensions. That is why each route needed a fresh check.

## Unifying law (exact, `d_price.py` and `verify.py` §1)
Let τ = Σ_a dim E_a be the number of unpaid port coordinates. Then mass = (73·112 − τ)·(8/9) + C per call. Break-even holds iff C < 64 + (8/9)(τ − 112).

For decoders whose terms are all 2-integral, Theorem A generalises (thakrarsagar, #1587): Σρ ≥ rank blockdiag Gram(E_a) = τ. So **any coverage-law change, whether uniform, non-uniform, paired or in any alphabet, needs dyadic (even-denominator) gain ≥ τ/9 + 35.6** over the integral floor τ.

| τ (labels) | break-even C | ×1.3 needs C ≤ | ×2 needs C ≤ | gain needed |
|---|---:|---:|---:|---:|
| 112 (lines, incumbent) | 63 | – | 28 | 49 (incumbent has 48) |
| 224 (planes/paired) | 163 | 152 | 127 | 61 |
| 336 (F8 lines) | 263 | 252 | 227 | 73 |

The bit side is the same to first order (111u ports, w = 888/9): gain needed ≈ 12.3u + 35.7.

The only known sub-unit mechanism is the E8 dyadic Gram G/2. It reaches 0.571 per coordinate. Every route below lives or dies on whether an analogous mechanism exists for a denser mask.

## The four routes, ranked

**1. (b′) F8-line ports, u = 3: open, but the evidence is negative.**
- F2^9 = F8^3, with dot = Tr Σx_i y_i (self-dual basis (3,5,7)). Both suppliers use 9-dimensional labels.
- An F8 line is a 3-dimensional unpaid subspace. Ports a and b are free iff h(a,b) = 0.
- Across 12 structures, the 112 ports fall on 58–62 F8 lines (72–80 port pairs share a line). There are 676–728 free pairs. α(Γ) = 40–43, proved optimal by CP-SAT (`b_graphs.py`).
- Floors: integral C ≥ 336; any decoder C ≥ α+1; PSD decoders C ≥ 3α ≈ 120. None of these closes the route.
- Numerics (uncertified, `lm_results.txt`): masked rank ≈ 45–46. Beating break-even therefore needs an average term cost ≤ 5.7. That is the loosest per-term budget of any route (planes need ≤ 4.2).
- In the natural cheap class, every term is supported inside one F8 plane: cost 3 if the plane is degenerate (9 such planes, one isotropic line Σx = 0) and 6 otherwise.
  - LM fails at every cost ≤ 390.
  - It first succeeds with all 64 nondegenerate planes plus 2×9 degenerate ones, at cost ≤ 438.
- F4 is impossible on F2^9, because the dimension is odd. On an even template, F4 labels are a special plane assignment, and PR2017 already proves that binary decoders there have gain 0 (C = 40 = naive).

**2. (d) Non-uniform upgrade: reduces to the full masked problem; the cheap classes are dead.**
- The bar is (8/9)(u−1) per upgraded port.
- **Decoupled designs** I_S ⊕ (G/2)[S^c] are always mask-valid, for any labels on S. They cost dim E_a = u per port, and u > (8/9)(u−1) for every u. So they always lose (`verify.py` §5).
- **Additive corrections** on top of G/2 fail. For an orthogonal 8-frame, ΔC ≥ 16 against a bar of 7 (our `frame-upgrade/`), and every upgraded plane port loses ≥ 24 of its 52 G-neighbours.
- What remains is a non-additive rebuild, which is the same open problem as route 1 or plane ports.

**3. (a) Paired families sharing a plane: dominated.**
- Sharing q_j as family i's second coordinate is exactly the plane assignment E_i = E_j = ⟨q_i, q_j⟩ (PR1825/PR2126). So it is not a new interface.
- On PR1825's matching, α(Γ_E) = 48 (exact). Generic planes reach 23. Masked rank ≈ 50–56 numerically, against ≈ 39 for a23. Pairing is therefore strictly worse than an unconstrained plane choice.
- **PR2126's C = 112 is a pricing artefact.** It charges a rank-2 coefficient matrix (I_2 on two ports) as a single term at support-pairing rank 2.
  - Under that rule, the line incumbent's G/2 would be one query over all 112 ports costing pair_rank = 9. That would give ceiling ×2.55 for free on today's word (`verify.py` §2).
  - In the rank-1 term model (one retained total per copy tree, #446 Thm C, which reproduces 64), the fit costs 4 per pair = 224 = naive. `ceiling.py` rejects that (no positive saving).
  - This is worth a review comment on PR2126 (not posted).

**4. (c) Output-side free coordinate: closed.**
- ctoth's rank invariant (#1587, comment 18914510 §2) lets only call crossings move the per-coordinate rank without payment. Each crossing moves it by at most the label dimension.
- So an "output-side" coordinate must be either a larger label (routes 1–3) or more calls. ctoth priced the more-calls option at −15.5% (quad template, m = 1 + 8N).
- No new event type is on the table. Its stop test is to name one.

## κ caps
- κ binds at min(complex, bit). Today the bit endpoint 0.0017334 caps any complex-only change at +1.1%, so every route must also run on the bit side. The bit side uses 9-dimensional labels, so F8 is available there too.
- Favourable caps on the complex side: u=2 at C=127 gives ×2.02 (κ ≤ 0.00341); u=3 at C=227 gives ×2.01; u=3 at C=252 gives ×1.3.
- Fully charged realised decoders, so far:
  - paired planes: 224 (fails);
  - F8: 336 naive (fails); plane-supported ≥ ~438 numerically (fails).
- No route is currently above the incumbent.

## Next experiment (ranked), with stop tests
1. **F8 dyadic lift.** Search for real vectors v_a ∈ R^k, k ≤ 46, with v_a·v_b = 0 on every non-free pair (an orthogonal representation of the F8 mask). Then N = Gram/|v|^2, and the terms are the coordinates. Price each coordinate by the F2 pairing rank of its support span. Mix line, degenerate-plane, plane and full supports in LM with a cost cap of 252.
   **Stop:** fail on 3 F8 structures × 10 restarts at cost ≤ 300 (well above the target). Then retire F8 and with it the whole alphabet route.
2. **Plane- or F8-specific bound ≥ τ − (τ/9 + 35.6).** This must charge cancellation between terms, since barrier B1 caps line-valid methods at 128 (planes) and 192 (F8).
   **Stop:** if no candidate inequality survives the a23 and F8-trial-6 instances in 1 session, stop.
3. Post the PR2126 reductio as a review comment (minor finding, so a comment rather than a PR).

Files: `common.py`, `a_pricing.py`, `b_graphs.py`, `f8.py`, `c_lm.py`, `e_planes.py`, `d_price.py`, `verify.py` (run with `$S/expO/run verify.py`, about 1 min).
