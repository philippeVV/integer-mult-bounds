# Selected stock envelope at native bootstrap stage 2 (conditional)

**A scoped, conditional source-bound derivation. It changes no κ, admits no supplier and adopts no outer proof.**
Card `selected-stock-envelope` (priority revision 4). The base is main `0fdd7d82`.

The stage-1 parent coefficient is already covered: PR2084 gives 8194, PR2120 gives 6748, and the
`findings/selected-stock-stage1` fork series (Discussion #1166, comments 18932239, 18932494 and 18934711)
closes the A₀ constant. Native stage 2 is the first stage whose stopped parent calls a **stopped** completed leaf
A₁ rather than the base leaf A₀. PR2183's map says later stages "need their own preceding completed ordinary
contract". This package covers that step for stage 2 only.

## Bound parameters (`ASSEMBLY.json` `bootstrap[1]`, recomputed from a₀)

- s = 67018507072324140916999892981492927675859083/4·10⁴⁶ ≈ 0.001675462676808103
- L₁ = 4,785,575,287,428
- L₂ = 1,704,765,459,304,159,860
- H₂ = ⌊L₂δ₂/2⌋−4 = 3,917,000,013, with exact maximality checked

The producer is PR589 BIT, m=64. It has 28,425 physical roles, maximum child 42 and row coefficient
18355 = 8194+9909+252. The hash chain is selected-result → ASSEMBLY → PR589 EXPECTED. A₀ has m=575, child 529,
28 bits, halving degree 9 and coefficient 252.

## Width inventory of the A₁ calls made by the stage-2 parent

The sources for this list are PR387 `PROOF.md:9`, PR185 reviewed `:13–15, :40` and `three-stage-cover-rows.tex:44, 50, 65`.

- **Calls inside the borrowed body**, at a fixed atom width w=⌈e^s⌉ for the whole tree:
  - weighted leaves (atom interchanges) at width w;
  - peeled atoms and bounded stopped leaves at width w;
  - fixed cover routing at "width O(w)", i.e. ≤ c_R·w.
- **Calls outside the borrowed body**, made after restoration: the high-part exchange at width k₂, and the remainder
  of fewer than w digits.

**Branch lemma (exact).** s·L₂ − L₁ ≥ 2,851,485,324,488,315. So for every e ≥ 2^L₂, every width-w call and the
exchange (k₂ ≥ 4096(w−1) ≥ w) run A₁'s **stopped** branch. Only routing calls of unstated width and the remainder can
reach A₀. Under the internal reading, stage 2 therefore needs exactly this: stage-1 fit at all widths ≥ 2^L₁ (the
PR2084/2120 stage-1 guard, which is conditional on their interfaces) and A₀ completed below that.

## Envelope

Write ν₁(x) for log_q of the reserve that an A₁ call of width x needs from its caller. All three readings satisfy
ν₁(x) ≤ ρ₁⌈x^s⌉(1+log₂x) for all x ≥ 1, with ρ₁ = 8194+63/10 = 82003/10:

- **Internal** (PR387:7, PR185:38): ν₁ = 0.
- **Retained** (A₀'s row is external): ν₁(x) ≤ (63/20)(1+log₂x).
- **Fully hoisted** (A₁ borrows nothing itself): ν₁(x) ≤ log_q S₁(x) + 2·(63/20)(1+log₂x). This uses PR2084's
  stage-1 bound, A₁'s own body calls at width ≤ c_R⌈x^s⌉ ≤ x, and the exchange/remainder < x.

`rows.tex:5,14` gives one reused inherited leaf reserve and W(w)^D. With PR2084's log_q W^D ≤ 8194·w(1+t),
where t = log₂e:

```
log_q S_2(e,w) <= 8194 w(1+t) + rho_1 ceil((c_R w)^s) (1 + log2(c_R w))
               <= (8194 + eps_2) w (1+t),
eps_2 = rho_1 * 2^(1 + s(gamma+1) - s(1-s)L_2),   gamma = log2 c_R,
```

The second line needs:

- ⌈y⌉ ≤ 2y;
- 2^{st} ≤ w ≤ 2^{st+1};
- γ+1 ≤ (1−s)L₂, so that log₂(c_R·w) ≤ t;
- γ+1 ≤ (1−s)L₁, so that c_R⌈x^s⌉ ≤ x inside A₁.

The checker confirms all of these side conditions exactly for every γ ≤ **γ_max = 4,777,557,234,645**, and confirms that γ_max+1 fails them. Values:

| γ | log₂ ε₂ ≤ | σ₂ |
| --- | --- | --- |
| ≤ 64 (e.g. c_R = 4096) | −2,851,485,324,488,301 | ≤ 8195 |
| γ_max | −2,851,477,319,869,468 | ≤ 8195 |

The parent stock coefficient is σ₂ = 8194 under the internal reading. Under the retained or fully hoisted readings it
is σ₂ ≤ 8194+ε₂. This replaces the earlier fully hoisted 2·8194+63/20, which assumed leaf width ≤ e. Without the
routing-width lemma (using only x ≤ e), the hoisted reading gives σ₂ ≤ 163943/10.

**Guard comparison.** ⌈log₂(8·8195)⌉ = 17 ≤ H₂ = 3,917,000,013. The fallback 163943/10 also passes. As a side check,
the conservative stage 2 meets the same conditions (L₂ = 1,704,765,476,910,090,189, H₂ = 3,917,000,033). The
headroom plays no part in the derivation. It is only the final comparison.

## Precise missing source lemma

**(W₂).** At atom width w, every fixed-cover routing call that the selected PR589 GL₆₄ cover issues to A₁ has width at
most c_R·w, with log₂c_R ≤ 4,777,557,234,645. Any radix-conversion factor at the same eligible q is included in c_R.

The sources give only "width-O(w)" (`three-stage-cover-bit.tex:208`, PR185:15). They were written for the historical
m=67 cover. A natural value is c_R ≤ 4097, the size of one G₁ block. The bound needs only the astronomically weaker
inequality above. W₂ is needed only for the retained and hoisted readings. The internal reading needs no c_R.

## Still conditional (not discharged here)

- PR2084's packet/compiler, fresh-field and restoration interfaces.
- A₁'s internal borrowing, the completed contract of PR387:7.
- H2 (time).
- Whether the m=64 cover realizes the word's operations.
- REM, which is not replayed.
- PR36 DAG regeneration.
- Stages 3–7 and the conservative chain beyond the side check.

The same argument applies at stage j≥3 with L_{j−1} and L_j, but this package does not check those stages. No full
cold replay, formal proof or independent human review is claimed.

## Reproduce

`python3 -I -B verify_stage2.py <repo-root>` prints `PASS_SELECTED_STOCK_STAGE2_CONDITIONAL_ENVELOPE`. It finishes
in about 1.5 s using the standard library and exact Fractions/ints. It checks 21 file:line citations and the hash
pins, and rejects all 10 mutation controls:

- the stage-1 cutoff used as stage 2;
- γ_max+1;
- ρ₁ one bit short;
- L₂−1;
- a false halving with child 46;
- 28426 physical roles;
- A₀ halving degree 8;
- a guard above H₂;
- a wrong citation line;
- importing 57604.

Credit: PR2084 and PR2087 (Adam Gall), PR2120 (Nick Petalas), PR2105 (Sricharan Vagvala), PR1956 and PR1965
(Josie Elliston), PR1990 (ToikY), PR589 (DreamingOfClouds), PR185 and PR387 (Rohan Arun) and the PR141 and PR36
lineage. Prepared with substantial Claude (Anthropic) assistance. Apache-2.0.
