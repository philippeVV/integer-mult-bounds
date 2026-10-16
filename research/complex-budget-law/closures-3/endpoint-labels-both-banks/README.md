# Endpoint labels on both banks (idea 8): closed in template form and, at every known decoder, in one-block form (×1.0)

Status: scoped closure, author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`.
Multipliers are over κ = 0.0017.

## Verdict

| Item | Value |
|---|---|
| Template form (the form expL5 priced at ×1.97) | **closed**. ctoth's rank invariant (#1587 c18914510 §2: U ≤ v(N − d) unpaid coordinates per replica, #1094 at equality) forbids any extra unpaid line coordinate at N = 9. The ×1.97 point is infeasible; the feasible realisation (more calls; ctoth's template) prices at ×0.845. Conditional on ctoth's hand proof, which was not re-derived here |
| Oracle (expNN's exact S_h-equivariant oracle on the #109 motif) | **strictly worse with any Y label**: all nonzero equivariant Y-source labels give s = Wm exactly and lose the motif's whole deficit (152 at h = 12; 16,560 at h = 40; 303,400 at h = 100; bit interface 6,680 / 141,700). The optimum is always "no Y label" |
| One-block form, every *dense* decoder | **closed** (Theorem W3): the both-bank threshold falls back from γ < 2/3 to the one-bank γ < 1/2. At every known decoder (γ = 4/7 on P112, 0.538 on P104, 0.6 on P120) the envelope is δ ≤ 0 |
| Exhaustive small instances | no instance with L < Σu (10 instances; data drops allowed) |
| Open | one named, unconstructed object: a decoder pair with ≥ 1/7 flat rows in two consecutive windows and cost < (2/3)V (§ Residue) |
| Multiplier | **×1.0** realised. Bit side: the oracle gives ×0 on the motif, so even a complex-only gain is capped at ×1.011 |

## Changed hypothesis

Y sources start at label 0 (03-motifs; Theorems W and B). Idea 8 gives each Y_b a source label ⟨q_b⟩ (q_b·q_b = 1, q_b ⟂ p_b) and
puts X_b's sink at q_b^⊥. The deficit identity D = 2(Σu − L) counts all source labels, so Σu goes from V to 2V.

## Theorem W2 (two banks; expNP's DM model, shears w1 = X→Y, w2 = Y→X, w3 = X→Y in one block F₂^h)

**L ≥ |D_X| + |D_Y| + C_p^{w1}(Q) + C_q^{w2}(Q) + C_p^{w3}(Q)**, with D_X, D_Y the ports whose data wires drop and
Q = P ∖ (D_X ∪ D_Y).
- Frames only climb on non-dropping wires, so X_b's gates have frames in [p_b, q_b^⊥] and Y_a's in [q_a, p_a^⊥].
- In w2 a drop-free path Y_a → X_b forces q_a ⟂ q_b, so the lossy part of the Y→X block is −I on the non-orthogonal pattern of
  {q_a}, which contains the diagonal. Scratch totals carry it, and each total loses ≥ rk G_q[R, S] (Theorem C).

Without data drops: L ≥ 2C_p + C_q, Σu = 2V, so δ ≤ 2(2 − 2γ_p − γ_q)/((2 + H/V)h). One-block threshold γ < 2/3 (one bank: 1/2).
Geometry does not raise γ_q: on P112, q_b := p_σ(b) with σ an orthogonal perfect matching satisfies every constraint
(`match_p112.py`), so γ_q = γ_p.

## Theorem W3 (one block, both banks, no data drops)

**Lemma I (inflation).** If T₁(a, b) ≠ 0 in the w1 window, a drop-free path X_b → Y_a exists, so frame(Y_a) ⊇ ⟨q_a⟩ +
span{p_b : b ∈ T₁(a)} from then on. Call a row *dense* if that span is p_a^⊥. For G/2 on P120, P112 and P104 every row is dense;
a row stays non-dense only if it keeps ≤ 32 of its free pairs (`flat.py`). By rank-8 rigidity
(`../../closures-2/p120-rank8-rigidity/`) every rank-8 P120 decoder is dense in every row.

**Lemma S (saturated sources).** Let B be the Y wires saturated at p_a^⊥ entering w2, and B₂ = {a ∈ B : σa ∈ B}. Then
L_{w2} ≥ |B₂|: a free entry (c, a) with a ∈ B forces q_c = p_a, i.e. c = σa; on a σ-transversal A the lossy block is I_{|A|};
a term touching A and σA gathers frame F and must drop ≥ 2; a term touching one block costs ≥ 1.

**Theorem W3.** L ≥ C_p(w1) + max(C_q(w2), |B₂(w1)|) + max(C_p(w3), |B₂^X(w2)|). If w1 is dense, L ≥ C₁ + V + C₃, so δ > 0
needs γ < 1/2. At γ = 4/7 on P112: L ≥ 64 + 112 + 64 = 240 > Σu = 224. At most one of two consecutive windows can use a dense
decoder without paying V.

**Lemma F (isolated up-set).** Frames whose perp contains no port label form an up-set 𝒰 that no data wire can share a gate
with; the rank μ held in 𝒰 rises by ≤ 1 per climb in and falls by ≤ 1 per drop out, and ends at 0. On P112 every Gram-term
centre lives in 𝒰. This is the provable form of "restoration drops per leaf": a naive leaf-per-port tree pays ≥ V − 8 = 104
exits per window. The literal per-leaf and per-port lemmas are **not** proved (a leaf staying inside q_b^⊥ pays nothing; the
data-drop gap of Remark W2.2 stays open for interacting ports).

## Residue (open)

δ > 0 with three compressing windows needs γ₁ + max(γ₂, 1 − 2f₁) + max(γ₃, 1 − 2f₂) < 2, where f_i is the fraction of flat
(non-dense) rows. At γ = 4/7 this needs f₁, f₂ > 1/7: at least 16 rows of both the w1 and w2 decoders must each zero ≥ 20 of
their 52 transfers, at cost < (2/3)V. None is known (rank 8 forces T = F; expNO's rank-9/10 fits stay dense; the structured
sparse screen needs rank ≥ 12). Even then, Lemma F's restoration exits and a bit twin come on top.

**Next step (cheap, decisive):** on P112, an LM fit at rank 9–11 with the flat zero pattern imposed. **Stop** if the flat
pattern's minimum rank is ≥ 12; then idea 8 is closed outright.

## Exhaustive small instances (`bfs.c`)

Exact bucketed Dijkstra over (frames, value matrix over GF(P)); moves are nested frame steps and transvections between registers
sharing a frame; data drops allowed. Scouting covered V = 1 (h = 2–4, up to 2 helpers, GF(2) and GF(3)) and V = 2 (h = 3, all
6 port configurations up to symmetry, GF(2), up to 9.9·10⁷ states): **none** with L ≤ Σu − 1. `bfs_small.sh` replays the two
smallest. Caveat: these toys have γ ≥ 1, so they test only the data-drop and restoration loopholes, not compression; a GF(p)
bound covers only schedules whose coefficients reduce mod p.

To replay a large instance (needs about 2.5 GB at log2cap 28):
```
gcc -O2 -o bfs bfs.c
python3 -I -B mkinst.py 3 2 0 '1,2;2,1' inst.txt && ./bfs inst.txt 3 28
```

## Pricing (`price.py`)

| setting | complex | bit | overall (κ = min) |
|---|---|---|---|
| template, expL5 (m′ = m + u, τ′ = τ) | ×1.97 | ×1.00 | ×1.011, and infeasible (rank invariant) |
| template, coverage law τ_unpaid = 224, C = 64 + 64 | ×2.00 | – | infeasible at N = 9 |
| template, Y window private (C_q = 112) | ×0.65 | – | ×1.0 (not used) |
| template, feasible (ctoth quad template) | ×0.845 | – | ×1.0 (not used) |
| one block, W2 envelope, P112, H = 0 / V / 8V | ×6.4 / ×4.3 / ×1.3 | needs its own bit word | superseded by W3: δ ≤ 0 at every known decoder |

## Files and reproduction

Run each from this directory. Runtimes under a 4G cap, one worker.

| Script | Checks | Runtime | Output |
|---|---|---|---|
| `uv run -q --no-project --python 3.12 --with numpy python -I -B eqv_motif_y.py 12` | oracle at h = 12: every nonzero Y label gives s = Wm; optimum = recipe | 8.5 s | `eqv_y_h12.out`, `PASS_EXPQT_EQV` |
| `uv run -q --no-project --python 3.12 --with numpy python -I -B eqv_motif_y.py 40 neg` | same in the bit interface at h = 40 (h = 100 runs the same way) | 20 s | `eqv_y_h40_neg.out`, `PASS_EXPQT_EQV` |
| `uv run -q --no-project --python 3.12 --with networkx python -I -B match_p112.py` | P112 orthogonal perfect matching (γ_q = γ_p) | 1.0 s | `match_p112.out`, `PASS_EXPQT_MATCH` |
| `python3 -I -B price.py` | W versus W2 envelopes, expL5 template formula, coverage law, supplier min | 0.1 s | `price.out`, `PASS_EXPQT_PRICE` |
| `uv run -q --no-project --python 3.12 --with networkx python -I -B flat.py` | Lemma I hypothesis: every G/2 row dense on P120/P112/P104; flat-row maximum 32 | 2.1 s | `flat.out`, `PASS_EXPQW_FLAT` |
| `uv run -q --no-project --python 3.12 --with networkx python -I -B replay_y.py` | both-bank replay engine (expNP's, with Y labels): Y0 L = 224 = Σu, Y1 L = 1008 | 1.7 s | `replay_y.out`, `PASS_EXPQW_REPLAY` |
| `bash bfs_small.sh` | builds `bfs.c` (`gcc -O2 -o bfs bfs.c`, into a temporary directory) and runs two V = 1 instances via `mkinst.py` | 2.0 s | `bfs_small.out`, `PASS_EXPQW_BFS` |
| `e8.py`, `e8map.json` | E8 line model (closures-2 copy) | — | — |

Sources: expQT and expQW. Prepared with Claude Code (Claude Opus 5.5) assistance. Not independently reviewed.
