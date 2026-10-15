# expSA — resynthesising the local word with a different gate multiset (helper lever)

Local scout, no κ claim. Model: selected #1094 k8 word (`../expBV/sel1094.json.gz`, `wordmodel.py` replay, which
reproduces local_profile/width_inventory exactly). Overlap check (`find_overlap.py`, 6 term sets): no prior work on
map-level / gate-multiset helper bounds beyond discussion #2278 (our own B1/B2).

## 1. Changed hypothesis
Keep the E8 decoder interface fixed (112 port lines, 8 centres at hyperplanes U_q, cut at target level (), the
post-cut correction α_t that each target must receive, the β target mixing) but let the gate list, gate multiset,
helper count and helper pins (start_extra / helper_final) change freely. B2 does not apply: it propagates frames
through *this* word's dependency DAG. The question is what any gate list implementing the interface must pay.

## 2. Moment bounds (repo `ceiling.py`, same call as propdag.py; `price.py`)
| scenario (selected word) | ceiling | κ if C scales (≈0.001715×) |
|---|---|---|
| any gate list with today's helper width inventory (Lemma 3 one-shot; width = dim end pin − dim start pin) | **×1.360** | ≤ 0.00233 |
| B2 (today's DAG, entry/exit free) | ×1.65 | — |
| zero helper mass, data registers charged at worst case (8 rank-1 climbs each, fully charged) | **×2.142** | ≤ 0.0037 |
| ideal: helper climbs removed, data as today | ×2.458 | ≤ 0.0042 |
| zero helper mass + one-shot data | ×3.47 | — |

Key point: subadditivity makes ×1.36 a cap for every gate list that keeps the inventory. The ×1.36 → ×2.46
range opens only if helpers become frame-constant, i.e. width 0 with start pin = end pin = working frame.
Frame-local touches between pinned ends do not get there. The bit side still binds at +1.1% (brief), so none of
this moves κ unless the bit side moves too.

## 3. Map-level lemmas (`deficit.py`, `classes.py`, `alpha.py`; every check passes)
- **L1, first touch (any gate list).** For each source s, the 52 sources its B-readers need are all ⊥ p_t, and
  the hyperplanes of its 52 readers meet in exactly ⟨p_s⟩ (112/112). Frames are monotone along every
  time-respecting path, so the first gate on s must be at frame ⟨p_s⟩. Before the cut, targets are frozen at (),
  and every other data register has its own port. So s's first partner is a helper, and the 112 sources need 112
  distinct helpers. The current word does exactly this (112 × (dim 1, helper)).
- **L2, island bijection (zero helper mass).** Frame-constant helpers at ⟨p_s⟩ can pass x_s only to registers
  sitting at () or ⟨p_s⟩. Before the cut, s must climb to reach its centres, to G_s ⊇ ⟨p_s,u⟩ with u ∈ K_s.
  After the climb at most **28** of its 52 readers satisfy u ⊥ p_t (for every u ∈ K_s). That leaves
  **2688 far (s,t) pairs** that must enter the target network through a target visiting the dim-1 frame ⟨p_s⟩.
  Each target visits at most one such frame. So zero helper mass forces an exact bijection targets↔sources
  (orthogonal), and every far pair must travel through target-to-target relays at frames of dim ≥ 2 that avoid u.
- **L3, A-side classes.** The 8 centre normals span 4 dimensions. The ports split into 14 syndrome classes of 8,
  with K_s = ∩_{q∋s} U_q of dim 6, and p_a ∈ K_b iff a and b are in the same class (784/784 vs 0/11648). A
  frame-constant carrier from a source's first read is class-pure, so with frame-constant helpers each centre can
  take first reads from one class only: ≤ 64 < 112. Class content must therefore leave K_σ through class-σ
  *source* registers that climb (data climbs, pinned width). This shapes the design but is not an obstruction.
- **Toy.** On F_2^3 with ports {1,2,4}, a hand-built schedule (`toy_check.py`, PASS) with one frame-constant island helper meets the
  analogous interface at zero helper cost. So L1–L2 alone do not rule out the ideal; the obstruction, if there is
  one, needs the 52-reader / 24-far-reader scale.

## 4. Answer and first obstruction
The fully charged reachable range is **×1.36 (proved cap for today's inventory) to ×2.14 (if a zero-helper-mass
word exists)**. That is 0%–100% of the "data-charged ideal" and is not narrowed further here. The first precise
obstruction: zero helper mass ⇒ the L2 island bijection plus relays carrying 2688 far pairs through targets
alone. Targets have ≤ 8 frame levels and each relay pollutes both frames.

## 5. Next experiment (ranked)
1. **CP-SAT on the L2 relay network** (frames as port spans, OR-reach relaxation, which is favorable): 112 targets,
   a bijection σ, relays at dim-2/3 frames. Decide whether all 2688 far pairs can be routed using targets plus
   *k* non-constant helpers, minimising Σ f(width). **Stop test:** if the dim ≤ 3 relaxation is infeasible with
   k ≤ 240 helpers of width ≥ 6, price the residual. If that is ≤ ×1.5, close the lever (κ×1.3 is unreachable
   here). If it is feasible, build the B phase literally and audit it with `audit_local`.
2. Re-run L1/L2 with β free (target mixing), since the α support, and so the 52-reader meet, depends on β.
