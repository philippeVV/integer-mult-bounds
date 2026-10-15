# Helper lever by gate-multiset resynthesis: the window is [×1.395 built, ×1.877 proved]

FINDING, no κ claim. This is a follow-up to `../helper-bounds/` (Lemma B2: with today's gate multiset
the helper lever is ≤ ×1.65). Here the **gate multiset is free**. The decoder interface of the selected
#1094 word is kept fixed: 112 port lines, the 8 centres at U_q, the cut at target level, the corrections
α_t and the β mixing. Pricing uses the repo's `ceiling.py` with data registers charged at their worst
case. Any complex-side gain is still capped by the bit endpoint, at about +1.1% over today.

The subdirectories are the steps of the chain, in order. Each has its own README and checks.

1. **`expSA`: interface and first cap.**
   - Structure: each centre is a ±1 sum of 56 sources. The ports fall into 14 syndrome classes of 8,
     and every target needs 52 orthogonal sources.
   - **Proved cap ×1.36 for any gate list that keeps today's helper widths.** Above that, a design needs
     frame-constant helpers.
   - L1: each source's first gate sits at ⟨p_s⟩ with a helper partner.
   - L2: a design with zero helper mass forces 2,688 "far" pairs to be relayed through targets.
2. **`expSF`: zero helper mass is impossible.**
   - Every target hyperplane needs a second register that arrives there, while the 8 centre hyperplanes
     use up 8 source slots. So at least 8 climbing helpers are forced.
   - With free centres, an exact all-ones split structure is forced, and a dim-6 count gives a
     contradiction: 2,912 > 2,688 (`cert.py`).
3. **`expSH`: per-level counting bound, ≥ 585 log-mass, i.e. ≤ ×2.031.**
   - Dim-2 lemma (`cert2.py`, tight locally): every port needs a charge of at least 3, so R₂ ≥ 336 and at
     least 112 helpers are alive at dim 2.
4. **`expSI`: an exact rectangle partition and a design.**
   - The 112 ports plus the 8 centre normals form 15 classes. Each class is the 8 odd vectors of a 4-dim
     space.
   - **All 6,272 demand pairs partition exactly into 392 rectangles odd(Q_A)×odd(Q_B)**, where Q_A ⊥ Q_B
     are 3-dim "quads" (`quads.py`).
   - α has rank 1 on the 336 cross-class rectangles and rank 4 on the 28 own-class ones.
   - A quad hub-tree design costs 2,893 log-mass (×1.709) in the OR-reach relaxation. Under the literal
     audit rules it costs **6,134 (×1.395)**.
     - R1: targets only absorb (no register may read a target).
     - R2: each edge carries one linear form.
     - The R1/R2 design is checked by `verify.py --r2` and `lincheck.py`, exact mod 2⁶¹−1.
5. **`expSK`: exact per-level MILP under R1.**
   - It covers all 8,282,955 subspaces of F₂⁹, classified by type (`prof.c`, checked by `chkprof.py`).
   - Result: **≥ 1,607 log-mass, so ≤ ×1.877**. If targets read only at their top, the bound is 2,368
     (×1.78).
   - The optimum suggests three changes:
     - targets absorb along their whole flag;
     - two line edges per port instead of three;
     - fatter middle nodes.

**Status.** The window is open in **[×1.395 built, ×1.877 proved]**. ×1.5 needs ≤ about 4,500 log-mass.
Neither a general closing certificate nor a ×1.5 design exists yet. An absorbing-target redesign is
running. A literal gate list would also have to pay for helper-junk cancellation, so ×1.395 is
optimistic.

**Run** (from each step directory, Python 3, standard library unless noted):
- `python3 expSF/cert.py`
- `python3 expSH/cert2.py`
- `python3 expSI/quads.py`
- `python3 expSI/verify.py --r2 <design>` and `lincheck.py`
- `expSK/chkprof.py`
- the MILP `r1ilp.py`, which needs highspy

Scripts that read the selected word use `../../helper-bounds/sel1094.json.gz`; generate it with
`../helper-bounds/run_all.sh`.

Prepared with substantial Claude Code assistance (sub-agents did the work; the main session re-ran the
key certificates from this directory). Not independently reviewed.
