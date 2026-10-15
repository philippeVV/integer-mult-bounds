# expM1: the framework's own ceiling

*Local scouting note. Nothing here was pushed or posted, and there is no new κ.
Base: main `0fdd7d82`, budget law `imb-work6/research/complex-budget-law`. The script uses only the standard library:
`python3 -I -B framework_ceiling.py` (about 5 s).*

## 1. Idea and changed hypothesis

Earlier ceilings fixed the local word (packing +0.87%), the template (levers.py) or one lever (budget law).
This note instead treats the whole recursive-supplier model family as a parameter space:
child size m, port label dimension u (unpaid coordinates per family), stages s, label dimension h,
decoder/centre cost C on V ports, local-word share φ and child ranks. It then maximises the exact
budget law over every choice the template identities allow.

**Changed hypothesis:** no particular word, decoder, lattice, stage count or port dimension is fixed.
Only the following structural identities are kept:

| # | identity kept | source |
|---|---|---|
| I1 | tilt identity: Ψ(a)=1 ⇔ E_π e^{aL} = 1 + D/M | budget law, Thm 1 |
| I2 | coverage: paid mass = (m−u)W + C·N; families per call w=(s−1)V/s ⇒ ν = D/W = u − (s/(s−1))·C/V | PR1012, #1587 |
| I3 | template: m = u + s(h−1) (selected 73 = 1 + 9·8) | template |
| I4 | **data floor:** each of the 2V data registers in a call climbs h−1 inside the stage, at rank ≤ h−1 (touch-separated). Share of paid mass φ_data = 2V(h−1)/(w(m−u)) = **2/(s−1)** (exactly 0.25 at s = 9, which matches the measured 1792/7232) | measured (budget law §3), upstream Lemma 5 table |
| I5 | other children have rank ≤ R, where m² > 2R² (outer halving guard) | ceiling.py |
| I6 | u ≤ h/2: a port needs orthogonal free partners, and pairwise-orthogonal nondegenerate u-spaces in the h-dimensional port space need 2u ≤ h (ctoth #1587: at most 4 pairwise-free planes in F2^9) | upstream §3 labels, #1587 |
| I7 | s ≥ 3. A swap needs three transvection stages, and φ_data ≤ 1 also forces s ≥ 3 | upstream Lemma (x,y)→(−y,x) |

## 2. The formula

Jensen applied to I1 gives exp(a·ℓ̄) ≤ 1 + D/M. Per family, D/W = ν and M/W = m − ν. Hence

$$
\boxed{\;a^* \le A=\frac{\ln\frac{m}{m-\nu}}{\bar\ell},\qquad \kappa\le\frac{A}{1+A},\qquad
\nu=u-\frac{s}{s-1}\frac{C}{V},\quad m=u+s(h-1),\quad
\bar\ell\ \ge\ \frac{2}{s-1}\ln\frac{m}{h-1}+\Bigl(1-\frac{2}{s-1}\Bigr)\ln\frac{m}{R}\;}
$$

Per unit of mass, the cost (r/m)^{-a} decreases in r. The lower bound on ℓ̄ therefore puts every class at its largest
admissible rank. The bound is tight to 0.06% on the selected profile: the exact root is 0.0016754627, Jensen gives 0.0016765,
and the selected word with ideal idle mass gives 0.00169005 (`ceiling.py`), reproduced in section A.
Writing ε = ν/m (net unpaid fraction per family) gives A ≈ ε/ℓ̄. This is the whole story:
**a is the unpaid fraction divided by the mean log child ratio.** Today ε = (5/14)/73 = 0.0049 and ℓ̄ = 2.93.

**Supremum over the whole family.** Take C = 0, u ≈ h/2 and s = 3. Then ε → 1/(2s+1) and ℓ̄ ≥ ln(s+½), so as h → ∞

$$A_{\sup}=\frac{\ln(7/6)}{\ln(7/2)}=0.12305,\qquad \kappa<0.1096\ \ (\times 64\ \text{over }0.0017).$$

A grid over s ∈ [3,40], h ∈ [3,200] and every u ≤ h/2 (section G) puts the maximum at s = 3 with small-h rounding effects:

| range of h | κ cap | × |
|---|---:|---:|
| h ≥ 9 | 0.118 (h = 10, u = 5) | 69 |
| any h ≥ 3 | 0.134 (degenerate h = 4, u = 2) | 79 |

Every one of these needs a **free decoder** and half-dimensional port labels.

With a decoder floor, the cap drops:

| decoder | h ≥ 9 | any h |
|---|---:|---:|
| γ = ½ (forced at u = h/2) | 0.036 | 0.042 |
| γ = 4/7 (E8 line ratio) | 0.027 | 0.032 |

## 3. Results (script sections; multipliers over 0.0017)

**Selected template (m=73, s=9, h=9, V=112, u=1)**

| scenario | κ cap | × |
|---|---:|---:|
| selected word, ideal idle (calibration) | 0.00169 | 1.0 |
| local word moved to rank 9 | 0.00288 | 1.7 |
| C = 0 only | 0.00473 | 2.8 |
| helpers deleted (data floor only), C = 64 | 0.00593 | 3.5 |
| u = 2 at line ratio, C = 128 | 0.00332 | 2.0 |
| **every lever free (C = 0, helpers gone, ideal idle), u = 1** | **0.0165** | **9.7** |
| coverage identity alone, ignoring the data floor | 0.0370 | 21.8 |

**Line ports (u = 1) over codes and stage counts (data-floor regime)**

- E8 P112: s = 9 gives 0.0059. The best s is 5 (m = 41), giving 0.0069 (×4.1).
- The open centre lever C(P112) = 56: 0.0094 at s = 4.
- The C ≥ n² values for BW16 and Leech are heuristic, not proven floors. They give 0.0149 and 0.0121, and only at s = 3.
- D4 and E6 have ν < 0, and E7 gives 0.0024.
- **With u = 1, no code or template exceeds about 0.015–0.035, even with C = 0** (0.0346 at C = 0, n = 8, s = 3).

**Growth directions (h → ∞)**

- **u fixed:** κ falls like 1/(s·h·ln s). For example, h = 129 gives 0.00038 and h = 1025 gives 0.00005.
  Larger codes help only through ν → u (at most ×2.8 at E8). E8 sits at the optimum of this trade.
- **u = k(h−1) (k ≤ ½):** κ saturates at a constant that does not depend on h. The values below are free decoder / Haemers floor γ = ½ / E8 line ratio γ = 4/7:

  | k | s | κ caps |
  |---|---|---|
  | ½ | 3 | 0.110 / 0.028 / 0.016 |
  | ½ | 5 | 0.085 / 0.033 / 0.025 |
  | ¼ | 3 | 0.064 / 0.016 / 0.009 |

- **No parameter direction makes κ grow without bound.**

**Decoder compression needed (section H)**

| target κ | condition on decoder ratio γ = C/(uV) | condition on port labels |
|---|---|---|
| 0.01 | γ ≤ 0.49 (k = 1/8) up to γ ≤ 0.75 (k = ½) | none extra |
| 0.1 | γ ≤ 0.06 | k ≥ 0.45 |

The κ = 0.1 conditions are not met by anything in the family:
- At u = h/2, a port's only free partner label is E_a^⊥. So a set of V/2 pairwise non-orthogonal ports carries an identity block, and the Haemers/min-rank argument gives **γ ≥ ½**.
- For k > ⅓, the orthogonality graph is triangle-free.
- No compressing u ≥ 2 decoder is known at all. The best plane decoder has γ₂ = 224/224 = 1, which gives ν < 0.

## 4. Answers

1. **Is κ O(1/log m)?** It is worse than that. With u fixed, a ≈ ν/(m·ℓ̄) = O(u/(m·ln(m/r))), which is Θ(1/m).
   The m = 73 denominator is the binding factor, not the logarithm. The largest κ with every lever free at u = 1 is about 0.0165 in the template (×9.7) and about 0.035 over all codes and stage counts (×20).
   With two or more unpaid coordinates and a free decoder, κ rises to 0.11–0.13 at most (×64–79). With any decoder obeying the half-dimension min-rank floor γ ≥ ½, it rises to about 0.036 (×21).
2. **Scaling with h and bigger codes.**
   - With fixed u, κ goes to 0 as h grows.
   - With u ∝ h, κ approaches a constant: at most 0.11, about 0.03 under the Haemers floor.
   - No direction is unbounded. The family's formal cap is κ < 0.134 (0.1096 as h → ∞), and the realistic envelope is κ ≲ 0.03.
3. **The single feature that must change.**
   - The net unpaid fraction ε = ν/m must become Θ(1) rather than about 1/m. In practice this means a terminal label (the "free" part of each route's endpoint frame) that is a constant fraction of the label space, with a decoder that still compresses at that dimension.
   - That alone saturates at κ ≈ 0.1. The data floor (I4) and the u ≤ (h−1)/2 orthogonality (I6) are built into the gather/side-injection/scatter label geometry.
   - So κ → O(1) near 1 needs ε → 1. Then the free pointwise-linear work would be doing almost the whole transform, which is the original problem again. Labelled subspace networks with nested residual pricing cannot do this.
4. **Can this framework reach κ ≈ 0.1 or 1?**
   - **κ ≈ 1: no.** κ = a/(1+a) < 1, and κ ≥ 0.5 needs ε ≥ 0.29, which is excluded by I3 + I6 (ε ≤ h/2 / (h/2 + 3(h−1)) ≤ 2/11 for h ≥ 4).
   - **κ ≈ 0.1: only at the formal edge.** It needs s = 3, port labels of dimension at least 0.45(h−1) and γ ≤ 0.06. At k = ½ these contradict the γ ≥ ½ floor, and nothing remotely like them is known.
   - **Realistic ceiling: κ ≈ 0.01–0.03** (×6–×20), and only through large-dimension port labels (u ≫ 1) with a γ ≈ ½ decoder.
   - **With line ports (u = 1): about 0.0069 using known decoders and data-floor packing (×4); 0.0165 even if every lever is free.**

## 5. Next step (ranked), with stop tests

1. **Prove a decoder compression floor γ_min(u, n) for u-dimensional port labels** (min-rank/Haemers via the orthogonality graph, then the label-pairing rank charge). This is now the only quantity between the realistic and formal ceilings.
   - **Stop test:** if γ_min(u, n) > (s−1)/s·(1 − m(1−e^{−0.0034·ℓ̄})/u) for every u ≥ 2 in the data-floor regime (section H's ×2 column, γ ≤ 0.71–0.87), the large-port route cannot even double κ. Then this whole family is capped near ×4 (u = 1), and work should move to a different outer recursion.
   - #1587 (ctoth) already shows line ports are at their cap (no second unpaid coordinate at nine calls), and that plane ports are exact on small instances but need C below about 163.6. Our formula is the same lever generalised to all u.
   - Cheap first instance: compute α of the orthogonality graph for the best known plane assignments (the α = 23 state `p292.json`). This gives an immediate Haemers floor γ₂ ≥ α/V to compare with 4/7.
2. **Check whether I3/I4 are forced** (m = u + s(h−1), data floor 2/(s−1)). For example, test whether a template with s = 4–5 stages exists for E8 lines. Section C gives ×4 at the data floor but only ×1.9 with today's helper word, so it is worth doing only together with a helper-free local word.
   - **Stop test:** if no s ≠ 9 E8 template satisfies the label nesting, drop it.
3. Note that **the bit side obeys the same law** and sits 3.5% above the complex side. Every number here is a joint cap.

## Scope

- These are planning ceilings under I1–I7, not lower bounds for multiplication.
- I4 and I6 come from the selected template and upstream label geometry. A network that prices frames differently (non-nested labels, shared children, a different outer recursion) is outside the family.
- Jensen and fractional counts make every number favourable. Realisation, hosting and precision are not charged.
- The BW16 and Leech C floors are heuristic.

Prepared with Claude Code (Claude Opus 5.5) assistance; not independently reviewed.
