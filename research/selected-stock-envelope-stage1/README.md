# Selected stock envelope, bootstrap stage 1 (conditional)

**Scoped source-bound derivation. No κ change and no supplier admission.**
This addresses card `selected-stock-envelope` (priority revision 4) for one
stage. #1956, #1965 and #1990 left the selected coefficient σ uncomputed.

## The stopped producer is the m=64 bit supplier

At a bootstrap stage the stopped producer is the selected **bit** supplier: PR589,
m=64, physical 28425, maxchild 42, halving depth 2. It is not the m=73 complex
supplier. The sources:
- `bootstrap-depth-audit/sources/pr387/refine.py:66` names the chain saving
  `degraded_bit_saving`.
- `PR386-BIT-GEOMETRY.txt:163-166` says the bootstrap cutoffs use the bit
  supplier's charged two-level coefficient (≤ 2^512).
- The selected assembler asserts that bound for the bit record.

m=73 enters only through the saving choice and the complex moment. For
completeness, the m=73 alternative (halving depth 1) also passes; see
`verify.py`.

## Envelope at stage 1

Width domain: e ≥ 2^L₁, t = log₂e, w = ⌈e^s⌉ ≥ 1, q > 2^80.

1. **Per-level stock.**
   - Source: `three-stage-cover-bit.tex:211-212`, `W(w) = W₀ q^{m²(w−1)}`.
   - The selected cover lies in GL₆₄, so `log_q W(w) < 4096w + 1`. This needs
     the per-level role multiplier N to be below q.
   - Overcounting N as physical × (1 + 32m²·paid) gives 52 bits, well below 2^80.
2. **Depth.** For n → r⌊n/64⌋ with r ≤ 42, 2·42² < 64², so two levels more
   than halve n. Hence `D(n) ≤ 2⌈log₂max(2,n)⌉` with no additive constant.
   This is proved, and also checked exhaustively for n < 2^20.
3. **Parent stock.** `D·log_q W ≤ 2(1+t)·4097w`, so **σ_new = 2(64²+1) = 8194**.
   - This matches PR589's recorded external row term 18355 = 8194 + 9909 + 252.
   - It is the selected-producer analogue of the historical 57604, and does not
     import 57604.
4. **Inner borrowing.**
   - Completed inner leaves borrow from their own chunks (PR185; PR387
     `PROOF.md:7`).
   - At most one old-leaf reserve R₀ is reused (`three-stage-cover-rows.tex:5-7`;
     upstream `04-swap.tex:354`).
5. **Old-leaf reserve for a₀ = 384599/10¹⁰.**
   - PR141 ties this exact saving to the old ordinary leaf row: 28 role bits,
     coefficient 252.
   - That gives `log₂R₀ ≤ 252⌈log₂max(2,e)⌉ + 28c`.
   - c is PR141's unquantified "fixed constant". The 252 is bound only to A₀'s
     own reserve, never used as σ.

**Result:**
- **σ₁ ≤ 8194** if A₀ needs no external reserve, which is how PR387 reads it.
- **σ₁ ≤ 8198 + c** under the three-stage-cover-rows reading.

The stage-1 guard `log₂(8σ₁) ≤ H₁ = 6,562,783` holds iff
c ≤ 2^6,562,780 − 8198. The threshold is checked exactly: it passes at the bound
and fails one above it. With c = 0, `log₂(8σ₁)` is about 16 bits.

## Missing lemma (the one coefficient that is still unbound)

The retained old ordinary supplier (a₀ = 384599/10¹⁰) needs a preceding-row
reserve with `log_q R₀(e) ≤ ρ(1 + log₂e)` for an explicit ρ. Equivalently, PR141's
additive depth constant c needs an explicit bound. Alternatively, a source proof
that A₀ needs no external reserve would close it.
- `three-stage-cover-rows.tex:5-7` says only "logarithmic".
- PR141 says "plus a fixed constant".
- PR387 `PROOF.md:7` says borrowing is internal.
These three readings have not been reconciled.

## Also open

- Applying `W(w)` to the selected m=64 cover is an inheritance reading. The
  source wrote it for an older m=67 word.
- The role-multiplier count is an overcount, not a source inventory.
- Stages 2–7 are not covered.

## Reproduce

```sh
python3 -I -B research/selected-stock-envelope-stage1/verify.py
```

Expected output: `"status": "PASS_STAGE1_CONDITIONAL_ENVELOPE"`. The script runs
in about 2 s with exact integers. It checks 17 cited source lines by substring,
the ASSEMBLY → PR589 hash pins and the H₁ recomputation. It executes none of the
sources.

## Open next questions

- Bound PR141's constant c, or prove internal borrowing for A₀.
- Repeat the derivation for stages 2–7.

Credit: Josie Elliston (#1965), ToikY (#1990), the #1956 authors, and the
PR141/PR185/PR386/PR387/PR589 authors.
Prepared by Philippe Van Velzen with substantial Claude (Anthropic) assistance.
A sub-agent derived the envelope, and the main session re-ran the check and
spot-checked the cited bindings. Apache-2.0.
