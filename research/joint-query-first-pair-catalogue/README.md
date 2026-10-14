# Joint-query controller: every catalogue first pair loses a query

**Scoped exact obstruction. No κ change.** This is for card
`joint-query-controller` (priority revision 4). It is the changed-first-pair
variant, which differs from #2074/PR2079 (an extra pre-F0 read) and from PR2086
(a changed F0 logger).

## Question

The selected witness `h = (48I − B₁*B₁)B₂*e` is lost by the first route pair
F0/B₁,₂: `B₁h = 0`, `B₂h = 2112e`. Does replacing that pair with another one,
`(F_a, B_{s,t})`, avoid the loss?

**Domain.** All 9 × 28 = 252 pairs with 0 ≤ a ≤ 8 and 1 ≤ s < t ≤ 8, i.e. the
9 F words and 28 B words of `research/joint-query` CONSERVED-QQ16. The following
are held fixed:
- the logical maps (1) and (3), with routes D_u = ρ_u(p) and D₀ = I;
- the A/T/U96 sets;
- the retained ordinary tail and the full observation set.

## Result: no pair avoids the loss

- **Each pair has its own lost witness.** For pair (a,s,t), take
  `h_{s,t,k} = (48I − B_s*B_s)B_t*e_k`.
  - `B_s h = 0` exactly, because `B_s B_s* = 48I` and the routes are involutions.
  - `B_t h ≠ 0`: its identity coefficient is `2112·I` for all 28 (s,t), it has
    4,225 group terms, and all 96 relative routes are distinct.
- **The loss.** An initial state with `X_U = h` and dyadic A/T compensation goes
  through `F_a` and then `B_{s,t}`, ending in the same state the reviewed witness
  reaches after B₁. So the reviewed 16-event tail applies verbatim, and every
  allowed observation is zero. A nonzero selected input is therefore
  indistinguishable from zero, also over ℚ.
- **General rule.** A pair (a,s,t) loses information whenever
  `ker B_s ⊄ ker B_a ∩ ker B_t`. This holds for all 252 pairs.
- **The selected vector itself.** `B_u h₁₂` is zero only for u = 1; for every
  other u it is nonzero (4,320 terms, or 4,225 for u = 2). Moving the first pair
  rescues that one vector but always loses another.
- **Overlap with PR2093.** PR2093's 3-port vector, which every B_s with s ≥ 1
  kills, already covers the a = 0 slice. This package also covers a ≥ 1.

## Reproduce

From the repository root, run `sh research/joint-query-first-pair-catalogue/reproduce.sh`.
It takes about a minute and needs only the standard library.

**What it runs:**
- `pair_catalogue.py` rebuilds R (8×112) and ρ_u from the pinned
  `community-frontier1042-audit/selected-replay.tar.gz` (sha256 `22c3f4dc…`).
  It then computes the group-ring witnesses for all 28 (s,t) and replays the tail
  from `QQ16-SELECTED-WITNESS.json`.
- `trajectory_check.py` runs the literal maps (3) then (1) for all 252 (a,s,t)
  with exact Fractions. It uses random port-diagonal permutation routes (m = 7)
  as a stand-in for the cover.
- Mutation controls: `ctl_47` (48 → 47) and `ctl_noA` (B_a cancellation dropped)
  must both fail.

**Expected output:**
- `PASS_QQ16_PAIR_CATALOGUE_WITNESS`
- `PASS_FIRST_PAIR_TRAJECTORY_ALL 252`
- both controls rejected

## Not covered

- Changed T records or ordinary responses (E_t, K_U).
- Pairs outside the 37-word catalogue.
- Paid reads from ordinary to selected registers, and affine helper state.
- Physical or native legality of the swapped frames.
- The tail is the reviewed restricted coefficient replay, not a native program.

To defeat the transported witness, a paid read must see `B_t h` on `ker B_s`
for the pair actually used; PR2093 shows that a single B_t read does not.
According to #2074, even an ideal rank-8 output leaves only +0.0437% ceiling
headroom.

## Open next questions

- Find the cheapest paid observation that separates `ker B_s` for some pair,
  and price it against the ceiling.

Credit: the `research/joint-query` authors; 1337isnot1337 (#2074, PR2079);
adamgall (PR2086, PR2093); Ross Richardson (#2074).
Prepared by Philippe Van Velzen with substantial Claude (Anthropic) assistance.
A sub-agent derived the result, and the main session re-ran every script and both
controls on main `0fdd7d82`. Apache-2.0.
