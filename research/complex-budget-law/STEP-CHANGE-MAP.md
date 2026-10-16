# Step-change map (2026-10-15): which structural routes could multiply κ, and which are closed

This is a synthesis of a day of scouting. It makes no κ claim. Every ceiling below is a favourable
planning bound, not a lower bound for multiplication. "Closed" means that within the stated model and
hypotheses the route cannot reach the stated multiple. The base is κ ≈ 0.0017; the claims on the
board are at 0.001717.

**How κ is built.** κ ≈ f / ln m, where f is the free fraction of each recursion level: about 0.72%
today. κ = 0.1 would need f ≈ 35%. The outer turns b into b/(1+b), and its fused form gives about b.
κ is the minimum over the complex and bit suppliers, and κ ≤ a_bit holds structurally, because the
complex supplier consumes bit calls through its change-of-basis adapters. So the **bit supplier is the
true ceiling** (κ ≤ 0.0017334 today).

## Closed or capped (fully charged)

| Route | Verdict | Where |
|---|---|---|
| Packing / owner fills | certified class ceiling; reached by #2745 | PR #2169 (7w/9 theorem), PR #2679 fill ceiling |
| Retired-lender aliases, pre-life hosting, frame-refit aliases | exhausted on the #2588/#2696 words (at most 6 hostings; refits other than 653 lose) | #2679, frontier census |
| Entrance raises, donor kernels, cone repair | none on the T74 words | #2453 comment; `levers/` |
| Helper lever, same gate multiset | ≤ ×1.59–1.65 (Lemma B2) | `helper-bounds/` |
| Helper lever, free gate multiset | window [×1.395 built, ×1.877 proved]; < ×1.5 once junk cleanup is charged | `helper-resynthesis/` |
| Larger local children | single-gate merges realise 0; larger frame spaces lose once charged | scout |
| Composition, nesting, bootstrap, coupled suppliers | ×1.00; the bootstrap is worth ~C² | scout |
| Recursive / tensored local words; paying a per-level overhead once | ×1.000: Ψ_{S⊗S′} = Ψ_SΨ_{S′}; cross-level sharing needs slot coincidence (Kronecker overlap; absent on #1094/#589) or frame coincidence (= flat tower). Ideal ×2.77 / ×9.65 unreachable; generous residual ≤ ×1.016 | `closures-2/zll-recursion/` |
| Flat tower (n > 4 exchange pairs per replica, m fixed) | ×1.00 realised (n ≤ 4: shared stage needs a rank-8 repair; 2n+1-call words have rank-1 blocks); ≤ ×1.05 for an unknown e ≥ 2 parallel word. The earlier ×1.32 row is withdrawn | `closures-2/flat-tower/` |
| Sibling calls sharing helper registers | capped at ×1.41 (k = 2, perfect rate-½ joint code; ×1.33 with junk); realised ×1.00; Lemmas P/S, reader cuts tight on `absorb10` | `closures-2/sibling-helpers/` |
| Other base objects (Leech, BW16, D_n, E6, E7, K12, …) | E8 is the unique Gram-type optimum; BW16 window [240, 549] open, ≤ ×2.8 | scout |
| Other coefficient algebras (F3, Z[ω], Z[i], Z[ζ_k]) | E8/F2 optimal ("why p = 2 wins": tight equiangular systems); #2017's C=40 has gain 0 | scout |
| Shared centre path between siblings | ×1.00 in today's recursion (36/36 chart pairs have pairwise-distinct routes) | scout; the fanout-recipe variant is still being tested |
| Removing the two-supplier min | structural (κ ≤ a_bit) | scout |
| Masked line decoders | every floor reduces to C(P120) ≥ 72; the lift family is proved ≥ 72 for |S| ≥ 2 | #1587 follow-ups 6–7 |
| Sparse-transfer ports (idea 4) and the corrector-rank lever, line ports h ≤ 12 | ×1.00: every rank-≤8 P120 decoder is Λ(G/2)Λ⁻¹ (T = F; DRUP-certified finite step); structured sparse supports with d̄ ≤ 30 need rank ≥ 12; half-integral decoders have T = F | `closures-2/p120-rank8-rigidity/` |
| ω·h floor (C ≥ ω·h) as a family screen | **false**: P112 C ≤ 64 < 72, P119 71, E7−1 48 < 49, F_7 50 < 56. Use per-object floors + monotone gain, γ(Q) ≥ 1 − (V_F − c_F)/|Q| | `closures-2/omega-h-floor-false/` |
| Field-line (polar-space) large ports, k ≈ 1/4–3/8 | ×1.00: best exact γ = 1 in all 8 families (frame decoders); minrank = V/ω law is numerical only | `closures-2/field-line-large-ports/` |
| TM-native / masked-lag pass programs | ×1.00: whole-volume pass model is κ-blind (≥ log_k(2^{u+1}−1) passes); XOR = route in every exact small table | `closures-2/tm-pass-programs/` |
| Coverage-law changes (planes, F8, paired, output side) | need even-denominator gain ≥ τ/9 + 35.6; lines get exactly 48 | #1587 follow-up 5 |
| Bit-side T74 transplant | ×1.000–1.017 (no disjoint lifetimes) | scout |
| Small cross-field units, h ≤ 4 | ≤ ×4.5 (odd-sink ratio lemma, x ≥ 2/3) | scout |
| Braverman–He F9 seed as a bit primitive (as written) | δ realised ≤ 0: helper restoration kills it; 3-integral and diagonal frames give δ ≤ 0 | scout |

**Framework ceiling** (`framework-ceiling/`): κ ≤ A/(1+A), with A = ln(m/(m−ν))/ℓ̄.
- **Line ports:** ≤ 0.0165 (×9.7) even with every lever free.
- **Large ports:** about 0.03 (×20) realistic, provided a compressing decoder exists. The formal edge is κ < 0.1096.
- **κ ≥ 0.5** is impossible in this family.

## Still open (live bets)

- **Small cross-field units at h = 5–6.** Certified x ≥ 0.4027 at h = 5 and x ≥ 4/15 at h = 6, so the
  ratio-bound ceilings (η = 1) are ×15 and ×19. No unit with x < 1 is known below E8. At h_c = 5 with complement sinks the
  pair-scored planning ceiling is ×3.03 ('best'; ×5.1 'cal') at the LP floor, and a unit only pays if
  x ≤ 0.7275. The integer certificate attempt ran out of memory and the CP-SAT unit searches are UNKNOWN,
  so the cell is open (`closures-2/small-units-h5/`). h = 6 is untested.
- **Characteristic-3 frames / Braverman–He import.** Local scouting (not yet packaged) found every tested
  signed-label family full rank under a generic-realisation screen (s ≥ rank Z_G for any char-0 gate
  assignment), i.e. δ ≤ 0 on those families; treated as paused, not closed.
- **Helper-free coding networks.** A characteristic-2 coding gap (Hermitian-twisted Tanner codes on
  PG(2,q) minus a matching; rigid and obstructed at q = 16 in local scouting) would be a bit primitive with
  δ ≫ 0.0067. Open obligations: an exact causal schedule and a frame family surviving the
  generic-realisation screen.
- **Large-port decoders.** No compressing decoder is known at label fraction ≥ 1/4. Proven floors:
  γ₂ ≥ 0.161 and γ₃ ≥ 0.357 (minrank/α bounds; they do not use the false ω·h floor). Field-line
  families are closed at γ = 1 (above); non-field u-dimensional labels are untested. An exact F8 screen
  found nothing below the naive cost (×1.0 realised).
- **Retained-array pricing convention.** Asked of maintainers in #1166; it decides the QQ16 and
  shared-centre routes.

Prepared with substantial Claude Code assistance (many sub-agents; the main session re-ran the key
checkers). Not independently reviewed.
