# Plane-port decoder: formal model (expG, working notes — not a publication)

Status: working formalisation for the lower-bound half of ctoth's #1587 plane-port question.
Sources: #446 Theorems A/C (research/five-cube-coded-copy-seed42/sources/PR446-PROOF.md),
#717/#968 decoder model (research/complex-frame-pair-decoders/PROOF.md), #1516 per-centre bound,
#1417 (P112 two-plane dimension >= 18), format-limits C.4 (degenerate labels dominated).

## Data
- Ports: the 112 binary ports q_a in F2^9 of the selected canonical word (gauge44 CANONICAL-WORD
  `ports`, inside research/community-frontier1094-audit/selected-replay.tar.gz). Form: ordinary dot
  product. All q_a odd weight (q.q = 1). P112 graph: a~b iff q_a.q_b = 0 (2912 undirected edges, 52-regular).
- Plane label: E_a = <q_a, r_a>, nondegenerate. Normalisation: every nondegenerate plane containing
  q_a has a unique basis {q_a, r_a} with r_a.q_a = 0, r_a.r_a = 1. So r_a ranges over the 128
  odd-weight vectors of q_a^perp (q_a itself is not in q_a^perp). Choice space 128^112.
  (Degenerate planes are excluded: format-limits Lemma Q says a degenerate-label unit is dominated
  by its quotient. Totally isotropic planes cannot contain q_a.)

## Contract (ctoth #1587)
source x_a: E_a -> FULL; target y_a: 0 -> E_a^perp; repair X_(q+r) = C_E^2. Each data incidence pays
h - 2 = 7; coverage becomes mass = 71W + C N. Break-even vs incumbent (deficit w - 64, w = 896/9):
2w - C > w - 64  <=>  C < w + 64 = 1472/9 = 163.55...  (first order, rest of word unchanged).

## What the local word must realise (scalar identity, #446 E5 / Theorem S)
After the frame-0 preconditioner the invocation acts as y <- y + J x with J = I (each target gets
exactly its own source). Split J = L + N:
- L = lossless (drop-free) part. By Theorem A / Cor. A1 (valid for any labels), a drop-free path
  x_a -> y_b needs E_a ⊆ frame(y_b) ⊆ E_b^perp, i.e. E_a ⊥ E_b. So L is supported on the graph
  Gamma_E = {(a,b): a != b, E_a ⊥ E_b} (never on the diagonal: E_a nondegenerate). Note
  Gamma_E ⊆ P112 graph because q_a ∈ E_a.
- N = drop-carried part = sum over copy terms k of t_k d_k^T (retained total d_k^T x gathered
  drop-free, delivered through one copy tree to readers with coefficients t_k).
Requirement ("masked decoder", the D1(a) class): N_aa = 1 for every port a, N_ab = 0 for every
a != b with E_a not ⊥ E_b; N_ab arbitrary on Gamma_E (fixed by the lossless part, priced at 0 here —
favourable).

## Cost (Theorem C, any copy tree; = #1516's per-centre bound)
rho_k = rank of the F2 pairing between P_k = sum_{a in supp d_k} E_a and Q_k = sum_{b in supp t_k} E_b
      = dim P_k - dim(P_k ∩ Q_k^perp) = rank_F2 of the block matrix [E_a^T E_b]_{a in R_k, b in A_k}.
C = sum_k rho_k. Incumbent line ports: E_a = <q_a>, P112 coordinate decoder, C = 64.

## Class covered by any lower bound proved here
Nondegenerate planes containing q_a; star/tree copy terms with one retained total each (separate
trees, #446 Thm C per tree); masked identity above; lossless part free. NOT covered: data-register
descents, shared drops across trees carrying several totals, isotropic/degenerate labels, changed
ports/h, nonzero-cost lossless network (helpers), target pre-mixing beyond the frame-0 preconditioner.

## Easy facts
- rho_k <= 9; rho_k >= 2 whenever supp d_k ∩ supp t_k != ∅ (contains E_a x E_a, rank 2).
- At most 4 pairwise-orthogonal planes in F2^9  =>  Gamma_E is K5-free (P112 has 8-cliques).
- Trivial decoder N = I: 112 terms of cost 2 = 224 (ctoth's naive number).
- Upper bound recipe: I = (G - A)/2 on P112 (G = Gram of E8 roots, rank 8; A = odd-pair part).
  Keep A on Gamma_E (lossless), re-decode A_lost = A off Gamma_E with extra terms.

## Results (expG, 2026-10-14) — see verify.py, breakeven.log, gfamily.py, natrank.py, minrank_numeric.py
- Exact: every port keeps at most 28 of its 52 P112 neighbours in Gamma_E, whatever the planes
  (|Gamma_E| <= 1568); edge-maximising local search reaches 1120 (unproved maximum).
- Rigorous floor (any planes in the class): drop cost-0 terms (they only touch Gamma_E); then each
  term costs >= 1, every diagonal needs a term with a in R∩A (cost >= 2), and N[S,S] = I_S on any
  independent set S of Gamma_E. So C >= rank(N) + 1 >= alpha(Gamma_E) + 1.
- alpha(Gamma_E) is NOT forced large: certified optima alpha = 56 (edge-max planes), 31 and 33
  (alpha-minimising search). So this floor is <= 32 for the best planes found: no kill.
- Natural decoder N = I + A_E/2 (incumbent E8 G kept on Gamma_E): rank 112 (mod 1000003, hence over Q)
  on all 6 instances; literal 'G-terms + lost-edge correction' costs 8*9 + rank(A_lost)=112 >= 184 > 163.6.
- Numerical evidence only (ALS floats): real minrank of Gamma_E ~ 56 for both an alpha=56 and an
  alpha=32 instance (d=45 fails, d=56 residual 6e-5, d=85 exact). A sub-break-even decoder at
  rank ~56 needs average term pairing rank <= 162/56 ~ 2.9.
- Break-even by term rank (ceiling.py unpaid=2, centre copies replaced by T terms of rank rho):
  max C = 162 (rho=1,2,3,6,9), 161 (7), 160 (4,8): break-even ~162 regardless of term rank.
