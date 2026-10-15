# expSI: building the hub-tree design for the #1094 interface (helper lever)

Local scout only. No κ claim, nothing pushed. This directory is a copy of expSH; expSA, expSF and expSH are unchanged.
The new files are listed at the end.

## 1. Changed hypothesis
expSH was bound-only. Here the interface is the same (112 port lines, 8 centres at U_q, α_t targets) and an explicit
register network is built for it. It is checked first under expSF's OR-reach rules, then under the first two literal
audit rules that those rules leave out.

## 2. Structure found (`quads.py`, all PASS)
- The 112 ports and the 8 centre normals form 15 classes. Each class is the 8 odd vectors of a 4-dim space.
- All 6,272 ordered demand pairs (5,824 (s,t) plus 448 (s,q)) are partitioned **exactly** into **392 rectangles of size
  4×4**: odd(Q_A) × odd(Q_B).
  - Q_A and Q_B are 3-dim "quads": halves of class spaces, with Q_A ⊥ Q_B.
  - All 14 quads of every class are used.
  - Each class pair splits into 2 blocks; the own class is K_{4,4}.
- So the cover the brief asked for exists with **no CP-SAT needed**. The hubs collapse to quads (dims 3 and 6), and each
  rectangle carries 16 pairs.

## 3. Relaxed design (expSF rules: OR-reach, targets may relay; `design.py`, `verify.py`)
- **Gather side** (per class): line → same-class pair node → quad.
  - Pairings are chosen by a Fano blocking set: e1 on 3 lines; e2 and e3 on 2 lines each; e1, e2, e3 collinear.
  - The minimum is 3 pair nodes per port, so every source's dim-2 frame lies in K_s automatically.
- **Middle:** one private register path per rectangle, Q_A → X → Y → Q_B^⊥.
- **Scatter side:** the mirror image (dual quad → dual pair → top or U_q).
- **Roles:** sources take pieces that start at their own line, targets take pieces that end at their own top (max-weight
  matching). The rest are helper pieces, priced by width.
- **verify.py PASS:**
  - flags are monotone;
  - each source starts at ⟨p_s⟩ with its dim-2 frame in K_s;
  - each target ends at p_t^⊥;
  - all 5,824 α pairs and all 448 centre pairs are delivered.
- **Cost:** 220 helpers, inventory {1:4, 2:48, 4:4, 5:28, 6:87, 7:49}, **2,893 log-mass, ceiling ×1.709** (price3).
  - So in the relaxation the lever is **≥ ×1.5**, and the stop test is passed with room to spare.
  - The relaxed window becomes [×1.709 built, ×2.031 LB].

## 4. Lifting: first obstructions
**R2 (first obstruction): `audit_local` forbids a non-target reading a target.**
- Targets can only absorb, so they cannot act as relays. The expSF relaxation lets them relay.
- Under R2, `design0.json` delivers only 3,419/5,824 α pairs and 304/448 centre pairs (`verify.py --r2`).

**R1 (linear capacity):** every register edge carries one linear form. Edge multiplicities were set from the exact
coefficient ranks (`alrank.py`):
- the 336 cross-class rectangles have α of rank 1, and each gather quad emits a single form;
- the **28 own-class rectangles have rank 4**;
- the centre rectangles (centre A-content) have rank 3–4.

**Repaired design (`coding.py`):**
- R1 multiplicities.
- Targets carry only their final edge D → p_t^⊥. Each target reads the registers at D.
- `lincheck.py`, exact mod 2^61−1: every emitted form lies in its node's in-span, and 120/120 tops and centres can form
  their exact demand.
- `verify.py --r2`: all PASS.
- **Cost:** 424 helpers, {3:10, 4:31, 5:100, 6:216, 7:67}, **6,134 log-mass, ×1.395** (7 seeds: 6,134–6,153).
- The interval lower bound for this edge multiset (`helperilp_coding.json`) is 5,977, so better splicing cannot reach ×1.5
  on this edge set.
- Most of the cost comes from R2: relay work moves from targets to helpers. R1 (rank-4 own-class rectangles and rank-3/4
  centre rectangles) adds about 144 edges at dims 3–5.

**Not done (task 3 literal lift):**
- No gate list was written, there was no `wordmodel.py` replay, and the decoder identity was not checked.
- The next obstruction is helper-junk cancellation. Every helper input is a free column, so a target must difference-read
  each contaminated arrival.
- Uncomputing in time needs common frames above p_t^⊥, which brings back target or helper frame growth. Doing it via
  cleanup at FULL needs FULL end pins (+1–2 width per helper) or restored helpers with buffers.
- Under these rules ×1.395 is therefore optimistic.

## 5. Answer
- **Relaxed (expSF rules):** a cover exists at 2,893 log-mass → **×1.709**. The relaxed lever is ≥ ×1.5.
- **First two literal rules (R2 targets-absorb-only, R1 linear capacity):** the best design found is 6,134 log-mass,
  **×1.395**. That is above the ≈4,500 stop line, i.e. **< ×1.5 in practice** for this design family. Its edge-set lower
  bound is 5,977.
- This is not a proof for all designs: the general LB is still expSH's 585, which ignores R1 and R2.

## 6. Next experiment (ranked)
1. **Lower bound with R2.** Redo expSH's per-level counting with targets as absorb-only:
   - only sources (112) plus helpers relay at dims 1–6;
   - add R1 rank cuts (the own-class blocks have rank 4).

   If the LB reaches ≈4,500, the lever is closed at ×1.5. Stop test: an LB below 3,000 leaves the question open.
2. **Cheaper R1/R2 design.**
   - Let targets carry B^⊥ → D → top when both tops of D are targets (target→target reads are allowed).
   - Route centre content through dedicated per-centre forms.

   Stop test: no design below 4,500 after these two changes → close.

New files: hubs.py (stats only; its 76 MB pickle was deleted), cls3.py, classes.json, quads.py, quads.json, design.py,
design0.json, verify.py, alrank.py, coding.py, coding0.json, coding4.json, lincheck.py, helperilp_coding.json.
