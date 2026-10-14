Task ID / proposed scope key: complex-new-helper-inventory / f3b0895a4f9f1cc62e75
  refined to: (1) k8-leaver K-block rank-4 cold ends; (2) E1->E2 entrance promotions
Chooser: seed 3c585d1f4e91703422bf417f1ab21f40, track complex-mechanisms, catalogue rev 2
Base commit: main ef5f39403610d231a0c27e8a6b4f08c10209bcd6
Predecessor / dependency heads: PR635 supplier a6895a69f8ad98ab45777b03e9e6f48888c39f6f
  (research/complex-571-cold14-st2, byte-identical to main's vendored copy);
  research/gauge44-leaver8 construction chain (selected lineage), leaver cell word 58e14a02...
Concrete object and question:
  (1) blocks {0,1,2,11,12,13} of the k8 leaver word: is rank(LK) > 0 achievable with the
      cells' exact scalar identity and 480 helpers kept? (cold44 currently pinned
      [0,0,0,4x8,0,0,0] -> 32 cold ends; the no-triangle word had 56)
  (2) 136 E1 entering helpers: which promotions to E2 are legal without induced climbs,
      and what is the exact per-call paid-profile change?
Scope status: refined. Reason: both questions change the helper endpoint inventory
  (width 9->4 cold ends; width 8->7 promotions) of a new word, i.e. a new provider record
  outside PR1034's frozen catalogue, as the card asks; they are the card's "narrow to one
  changed boundary or helper type with explicit counts" variation.
Existing restriction being changed: #1005 chain's pinned pairing ranks; BlockInvert's
  promotion recipe limited to E4->E5.
Method: exact GF(2)/rational recomputation of L, K, rank(LK) per block; declared variation
  class (completion priorities, source-copy exchanges, carrier placements, cell
  relabellings); promotion legality generalised from build_promotion.py; native
  independent audit_local.py on any candidate word.
Nearest overlapping work: #1529 (partial leaver cells beyond the eight; pairing collapse
  observation), #1148 (E4/E5 promotions, E1 removal), #1138/#1126 (leaver relabelling
  atlas / 21-ADD minimum), #1240/#1272 (packing ceiling, rev-3 cards).
Optimistic full cost budget (frontier dual LP, relaxation): (1) +0.079%/end, <= +1.9% for 24;
  (2) ~+0.023%/clean promotion, <= ~+1.8% for 96. Required observable improvement: dRC < 0
  at the frontier dual and a profile change confirmed by the native audit.
Success evidence: audited word with rank(LK)>0 in a K-block / with a clean promotion.
Falsification test: rank(LK)=0 forced over the declared class / every promotion induces
  >= 1 climb.
Time/compute cap: ~3 h wall clock for both, sub-agents in parallel.
Predecessor authors: icekylinx (#645), PartialK (no-triangle), LeaverSAT/LeaverSAT2 (#1005
  lineage, jeremykothe2), BlockInvert (promotions), NewBArch, G44Pack3, higherbeing (#1529),
  Ali Romman (#1034 catalogue bound), Douglas Colkitt (#1240 ceiling).
Intended artifact: scoped construction or exact obstruction per question; PR only if it
  changes something maintainers must decide, otherwise a FINDING/HANDOFF post.
Claim: https://github.com/CrocSwap/integer-mult-bounds/discussions/1166#discussioncomment-18913759
Outcome: both questions answered negatively within their declared classes (see README).
  (1) dim L = 0 in every K-block for any circuit on a cell's own registers (pigeonhole lemma);
  (2) 104/136 E1 helpers have no legal E2; the 32 legal promotions each cost dRC = +0.0023566.
Stop rule used: exact obstruction or exact negative price; global repacking not run.
