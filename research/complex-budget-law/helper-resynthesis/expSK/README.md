# expSK: per-level helper lower bound under the real audit rules (helper lever, #1094 interface)

Local scout only. No κ claim, nothing pushed. This directory is a copy of expSI; expSA, expSF, expSH and expSI are
unchanged. The new files are listed at the end.

## 1. Changed hypothesis
expSH's per-level bound (585 log-mass, ×2.031) let all 224 data registers relay. Here it is redone under the two rules
expSI identified:
- **R1 (absorb-only targets):** a non-target may never read a target. So only the 112 sources and the helpers carry edges.
  A target t can only *absorb*: it reads non-target registers at the nodes of its own flag inside p_t^⊥.
- **R2 (one form per edge):** each register edge carries a single linear form.

The bound is in edge form, which is sharper than expSH's register form:
- E_d is the number of non-target edges at step d → d+1.
- The constraint is E_d ≤ 112 + (helpers covering step d). Sources are free, since they are charged worst case.
- A helper covers the step interval a..b and costs f(b−a+1), with f(w) = w ln(73/w).

## 2. Model (`r1ilp.py`, exact MILP with HiGHS)
1. **Node profiles (`prof.c`).** All 8,282,955 subspaces X of F_2^9 at dims 1–8 are classified by
   (P = #ports in X, C1/C2 = #ports covered by one/two hyperplanes of X, QV/Qp = #V/#ports in X^⊥, R1/R2 = #V covered by
   one/two hyperplanes of X^⊥), where V = ports + 8 normals. This gives 8–104 types per dim.
2. **Node capacity.** A node with k non-target in-edges, j out-edges and z absorbing targets carries
   ≤ c(k)·(r(j)+z) pairs, where c(k) = C1, C2, C2+C1, … (capped at P) and r(j) likewise from R1, R2 (capped at QV).
   z ≤ Qp, and r+z ≤ QV. The number of nodes of each type is capped by its subspace count.
3. **Absorption bookkeeping.** Each target sits at one node per dim, so there are ≤ 112 absorbers per dim. Cumulative
   absorption by dim d is ≤ 112·Mtop(d), with Mtop = 1, 2, 4, 8, 12, 20, 32, 52.
   - Pairs not yet absorbed pass a non-target node at every dim.
   - Pairs absorbed above dim d must leave their dim-d node through out-edges.
   - The 448 centre pairs (s,q) are never absorbed by targets.
4. **Arrival cuts (R2 / content).** Every top needs ≥ 1 non-target arrival (52 > 32). Every U_q needs ≥ 2 arrivals
   (56 > C1 = 32). Together, E_7 ≥ 128.
5. **What R2 adds, and what it does not.** The only R2 cuts that hold for *every* network are the per-node span and
   per-top arrival cuts above, and they are not binding.
   - The per-rectangle rank cuts (rank 4 for own-class rectangles, 3–4 for centre rectangles) are valid only for designs
     that route each rectangle through a private path, as expSI does.
   - A general network can feed a target the same sources along several paths, so these cuts do not hold in general and
     are **not** in the bound.
   - Their likely size: expSI measured about 144 extra edges at dims 3–5, so ≲ 340–620 log-mass. This is an estimate,
     not part of the proof.

Checks: `chkprof.py` (stdlib) **PASS**:
- the subspace counts equal the Gaussian binomials;
- dims 1 and 8 recomputed exhaustively in Python match;
- 150/150 random subspaces (dims 2–5) and 40/40 (dims 6–7) land on listed types;
- Mtop is as stated.

HiGHS reports optimal status on every run.

## 3. Results
| variant | min helper log-mass | helper intervals (steps) | ceiling (price3) |
|---|---|---|---|
| expSH, targets relay (reference) | 585 | 123×w1 + 8×w2 | ×2.031 |
| **R1, targets absorb at any dim (the bound)** | **1,607** | [1,3]×37, [1,6]×55, [1,7]×16, [2,3]×23 | **×1.877** |
| R1 without the absorption bookkeeping | 1,541 | — | — |
| R1, absorption only at dims ≥ 3 / 5 / 7 | 1,668 / 1,798 / 2,206 | — | — |
| R1, absorption only at the top (expSI coding.py regime) | 2,368 | [1,7]×112, [2,6]×21, [2,7]×16, 2×w1 | ×1.778 |
| expSI coding design (R1 + R2, built) | 6,134 | 424 helpers | ×1.395 |

Non-target edges per step (d = 1…7):

| | step 1 | step 2 | step 3 | step 4 | step 5 | step 6 | step 7 |
|---|---|---|---|---|---|---|---|
| LB | 220 | 243 | 243 | 162 | 183 | 183 | 128 |
| expSI coding design | 336 | 496 | 536 | 536 | 536 | 504 | 248 |

## 4. Answer
- **The lower bound under R1 (plus the network-valid R2 cuts) is 1,607 log-mass, i.e. the ceiling is ≤ ×1.877.**
- This is far below the ≈ 4,500 needed to close the lever at ×1.5, and below the 3,000 stop line. **The helper lever
  stays open:** the proved window is [×1.395 built (expSI), ×1.877].
- Even the strictest absorption reading (targets read only at their top) gives just 2,368 (×1.778). Adding the full
  estimated rank cost (≈ +600) still leaves it below 3,000.
- The gap to the best design is ×3.8 in log-mass (×2.6 in the top-only regime). Per step it is concentrated in the
  middle steps 3–6 (≈ 536 vs 160–245 edges) and in the line out-degree (3 vs 2 per port).
- No certificate is written: the outcome is "open", not "closed".

**Cheapest structural changes suggested by the MILP optimum:**
1. **Targets as multi-level absorbers.**
   - In the optimum, targets absorb at most dims from 3 to 7. All 112 absorb at dims 3–4 (2–4 ports each); 67–90 absorb
     at dims 5–7 (6, 15 and 20 ports).
   - R1 allows this, because targets only read. The coding.py design lets targets read only at the last edge.
   - Allowing it lowers the bound from 2,368 to 1,607 (≈ 760 log-mass, about 1/3). The middle edges drop most, since the
     scatter tree largely disappears.
2. **Two line edges per port instead of three.** 52 + 4 demands = 56 = R2 of a port line, so two out-edges suffice in
   the LP. The design pays 3 (Fano pairing, 336 vs 224 step-1 edges).
3. **Fatter middle nodes than the 4×4 rectangles.**
   - The LP uses nodes carrying 24–36 pairs per edge: dim-4 nodes with 6 ports and reach 6–8, dim-5 class-space nodes
     (8 ports → 4 normals) for the centre pairs, and dim-6 nodes with 15 ports and reach 3.
   - The design's rectangles carry 16 pairs per edge.
   - Under R2 this only pays where the α block on (content × reach) has low rank. That is the open question.

## 5. Next experiment (ranked)
1. **Absorbing-target redesign of expSI coding.py.**
   - Give each target a full flag inside p_t^⊥ through one gather quad Q ⊥ p_t (dim 3), its class hub, and the dim-6/7 pair
     levels. The target reads the non-target registers it meets. Keep the gather side and the 392-rectangle middle, but drop
     every scatter edge whose pairs some target already absorbed.
   - Verify with `verify.py --r2` and `lincheck.py`.
   - **Stop test:** if no seed gets below 4,500, the practical lever stays < ×1.5. If one does, lift it literally
     (gate list, `wordmodel.py` replay, helper-junk cancellation).
2. **Rank-aware LB (to close the gap from below).**
   - Replace node capacity c·r by the maximum nnz of an α-block with rank ≤ k, in the no-cancellation regime where each
     pair is served once. This needs per-node α data, not just types.
   - Add a charge for cancellation: overlapping supports cost OR capacity.
   - Only a ×2.8 jump would close the lever, so this is lower priority than item 1.

## 6. Discussion-post draft
> **Helper lever on the #1094 interface: lower bound under the real audit rules.**
> Headline structure (expSI): the 112 ports plus the 8 centre normals form 15 classes, each the 8 odd vectors of a 4-dim
> space. All 6,272 demand pairs (5,824 α pairs plus 448 centre pairs) partition exactly into **392 rectangles
> odd(Q_A) × odd(Q_B)**, where Q_A ⊥ Q_B are 3-dim halves ("quads") of class spaces. α has rank 1 on the 336
> cross-class rectangles and rank 4 on the 28 own-class rectangles. A quad hub-tree design built on this cover costs
> 2,893 helper log-mass (×1.709) in the OR-reach relaxation. Under the literal rules it costs 6,134 (×1.395): targets
> are absorb-only, and each edge carries one form.
> New: redoing the per-level lower bound with absorb-only targets (exact MILP over all 8.3 M frame types, edge form,
> absorption bookkeeping) gives **≥ 1,607 helper log-mass, i.e. ceiling ≤ ×1.877**. That was ×2.031 before. If targets
> may read only at their top, the bound is ≥ 2,368 (×1.78).
> So the lever is **not closed at ×1.5**: the window is [×1.395, ×1.877]. The bound's optimum points at the gap: targets
> absorbing along their own flags (quads, then class hubs), two line edges per port, and middle edges carrying ~30 pairs
> instead of 16. The next step is an absorbing-target redesign, with stop line 4,500. Even a full win here still needs
> the bit side to move (bit cap +1.1%).

## Files (new)
- `prof.c`, `prof1.txt`…`prof8.txt`: node profile types.
- `vecs.txt`
- `chkprof.py`: stdlib cross-check.
- `r1ilp.py`: MILP. Usage: `r1ilp.py R1 [amin=k] [nosplit] [-v]`; outputs `r1ilp_R1_*.json`.

All jobs were run capped, with at most 1 worker. No background processes remain.
