# expSF: is a zero-helper-mass gate list possible on the #1094 interface? (helper lever, x1.36 vs x2.14)

Local only, no κ claim. Inputs are copied from expSA (`instA.json`, `alpha.json`, `cen.json`, word `../expBV/sel1094.json.gz`).

## Relaxation used (strictly more favourable than the brief's CP-SAT model)
Gates join frames: both registers climb to the gate frame, which must lie inside both end pins. Framed x-content held by
a register has its ports inside that register's frame. Cancellation and time order are ignored (OR-reach). Every
helper may be **frame-constant at any pin, at zero mass**; such helpers act as hubs. Sources and targets act as relays
at all dimensions, and data climbs are free because they are charged worst case. Each chain can be refined to a full
flag at no loss. Under these rules the post-cut network is time-free. The content at a frame value F is the union, over
the registers that ever sit exactly at F, of what each one held just before it got there (cluster lemma).

## Results (all checks pass: `cert.py`, plus the exhaustive `exact56.py`, `max5in7.py`, `betafree.py`)
1. **Zero helper mass is infeasible. Certificate:**
   - α_t is exactly the 52 ports of p_t^⊥ and spans it.
   - Any 7-dim subspace of p_t^⊥ holds ≤ 32 ports.
   - So every target hyperplane (top) needs at least one contributing arrival besides t itself: a source or a
     non-constant helper whose frame passes through exactly p_t^⊥.
   - Each register has at most one 8-dim frame value.
   - None of the 8 centre normals is a port, so no U_q is a top. Each U_q needs its own non-constant register at
     frame U_q (a source or a climbing helper).
   - So 112 tops need ≥ 112 arrivals, but only ≤ 104 sources are left: **k ≥ 8 non-constant helpers, each of width
     ≥ 1.**
2. **It is not just the centres.** Even if the centres cost nothing, the 224 free registers alone cannot do it. The
   top count is then exact: each top gets t plus one source. This forces:
   - Each register's dim-7 level is a 26-port pair level ⟨p_t,p_u⟩^⊥. G26 is 48-regular; at each top the two levels
     split its 52 ports disjointly, and p_b = p_a + j + p_t with j = all-ones.
   - Each such level holds exactly 2 registers, so Σ_r c(r,6) ≥ 112·26 = 2912.
   - Every dim-6 level carries ≤ 12. A lone register is limited by its 5-dim level, which holds ≤ 10 ports inside a
     G26 level. A shared level has ≥ 4 ports in L6^⊥, so it holds ≤ 12.
   - 224·12 = 2688 < 2912, a contradiction.
3. **β free (L1/L2 redo):** for every target and every β, f_t − Σβ c_t' cannot be confined to a 7-dim subspace of
   p_t^⊥ (exact, mod p, `betafree.py`: 0/112). So the top double-arrival obstruction does not depend on β. L1's meet
   was not recomputed for free β; the obstruction does not need it.
4. **Construction search:** annealing on the hub model, using bijective σ/τ matchings and moves that snap chains onto
   levels already in use (`hubsa.py`, `hubsa2.py`), reaches only ≤ 279 of 5824 pairs. The design is very rigid:
   exact disjoint splits at dims 8, 7 and 6. Local search does not find it, and I have no construction. **Best feasible
   point known = today's word.**

## Pricing (`price2.py`, same ceiling call as expSA; data charged worst case; one-shot helpers)
| helper inventory | ceiling |
|---|---|
| none (unattainable, item 1) | ×2.142 |
| 8 × width 1 (**proved lower bound**) | ×2.136 |
| 112 × w1 / 112 × w2 / 112 × w8 | ×2.050 / ×1.995 / ×1.833 |
| 240 × w6 / 240 × w9 | ×1.631 / ×1.551 |
| today's inventory (one-shot) | ×1.258 |
| ×1.5 break-even | ≈ 4,500 helper log-mass (≈ 1,069 width-1 helpers; just above 240 × w9) |

## Answer
A zero-mass gate list does not exist, even in the most favourable relaxation. The proved minimum is tiny, though:
8 width-1 helpers, about 34 log-mass. That moves the fully charged ceiling only from ×2.142 to ×2.136. **The stop test
is not met.** Closing the lever at ×1.5 needs a lower bound of about 4,500 helper log-mass, and counting at the top
levels gives about 34 (the dim-6 deficit adds about 224 content units, roughly 19 width-1 helpers' worth if the
missing capacity were made up at dim 6, ×2.126). The lever stays open in [×1.36 proved cap, ×2.136]. There is no B
phase to audit literally, because no relaxed design was found.

**First exact obstruction.** Every target hyperplane needs a second arrival (t's own chain covers ≤ 32 of 52), and
the 8 centre hyperplanes, none of which is a top, use up 8 of the 112 source dim-8 slots.

## Next experiment (ranked)
1. **Exact minimum helper count/mass on the hub model, top-down.** Encode dims 8→5 as CP-SAT/ILP with clusters as
   variables:
   - tops: t + arrivals;
   - dim-7 levels: the G26 levels plus unpaired ones;
   - dim-6/5 levels: the hyperplanes of the dim-7 levels that are used (127 per level).

   Add content-count constraints c(L) ≤ min(ports(L), Σ pre) and helpers as extra cluster members costing f(width).
   Every bound in items 1–2 is a cut of this LP/ILP. **Stop test:** if the LP bound is ≥ 4,500 log-mass, close the
   lever at ×1.5. If it stays below about 500, the lever can only be settled by an explicit design: build it top-down
   (exact j-splits) and replay it with `wordmodel.py`.
2. Look for a symmetric construction. The ports are odd vectors and every 2-cover is an exact j-split, which suggests
   a recursive design that uses the symmetry group of the 112 ports. Then add the minimum number of helpers at the
   dim-6 deficit.
