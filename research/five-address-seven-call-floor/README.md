# Five-address dyadic model: the exact call minimum is 8

**Scoped exact obstruction. No κ change and no supplier change.**

This package covers card `bit-new-register-model` (track bit-mechanisms, priority revision 4, base `0fdd7d82`). The card came from the second suggestion of chooser draw `84c4cc9f…`. The first suggestions of four draws all landed on cards already covered in this session; see `draw*.txt` and `TASK-CARD.md`. It narrows the card to the open question left by PR #2112 (Discussion #2091): for that PR's five-address rational frame alphabet, is the minimum number of calls 7 or 8?

## Model

The setting is PR #2112's model:
- Coefficients live in F₂[Λ] with Λ = ZA + ZP + ZR.
- The frames are A = E₁₁, B = E₂₂, P = ½J, Q = A + B − P and R = (1,½)ᵀ(1,½).
- There are four dirty roles with free GL₄(F₂) mixing.
- The target is diag([A+P], [A+Q], [B+P], [B+Q]).

## Result: the minimum is 8, so the direct word is optimal

1. **First order.** The R calls must cancel at first order (S_R = 0), so the number of R calls m_R is even.
   - If m_R = 0, the word lives in the old lattice, whose first-order floor is 8.
   - If m_R ≥ 4, the remaining labels need at least 4 more calls, so the total is at least 8.
2. **The case m_R = 2.**
   - The homomorphism Λ → Z² that sends R to 0 induces a ring map on the group algebras, and every R call becomes I under it.
   - So a 7-call word would give a 5-call word over {A, P, Q, B} realising diag([(1,1)], [(4,3)], [(2,5)], [(5,7)]).
   - First order leaves only four profiles: (P2,Q2) with 6 calls, and (P2,Q2,B1), (A1,P1,Q3), (A1,P3,Q1) with 7 calls.
3. **Exhaustive exact search over F₂[Z²].**
   - 169,920 distinct 7-call residual words: 0 hits.
   - 216 distinct 6-call words: 0 hits. This also confirms #2112's six-call exclusion.

## Checks

| Script | What it does | Expected output |
|---|---|---|
| `seven.py` | Search, plus the direct-8 positive control and the drop-one negative control. | `controls ok`, `seven_hits: 0` |
| `check2.py` | Independent implementation: projectors built as u vᵀ, different enumeration, subset-expansion product. | Same counts, 0 hits |

Each run takes 1–3 minutes, needs only the standard library, and is invoked as `python3 -I`.

**Limits:**
- Only paid calls are counted, in a favourable formal model. Frame setup and additions are not priced.
- Not covered:
  - a fifth frame that survives every target-preserving projection;
  - targets that need R at first order.

Credit: Cinna Davis (PR #2112), 1337isnot1337 (#2091), and the #1974 and #2027 authors.

Prepared by Philippe Van Velzen with substantial Claude (Anthropic) assistance. A sub-agent ran the search. The main session re-ran both implementations, which agree. Apache-2.0.
