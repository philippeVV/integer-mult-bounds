# expSH: minimum helper mass in the expSF OR-reach relaxation (helper lever)

Local scout only: no κ claim, nothing pushed. This directory is a copy of expSF. The expSF and expSA files are unchanged;
new files are listed at the end. The model is expSF's hub relaxation: register chains are refined to full flags, frame-constant
helpers are free hubs, timing and cancellation are ignored, a source's L2 lies in K_s, and the 8 centres need registers at U_q.
Content then flows along the union DAG of register edges.

## 1. Changed hypothesis
expSF bounded only the top levels (dims 8, 7, 6). Here every level d = 1..8 is bounded by a counting argument:
- Each of the 5,824 pairs (s,t) passes exactly one node X per dimension.
- The pairs through X are at most content × reach, with content ⊆ ports(X) and reach ⊆ ports(X^⊥).
- The number of registers at X is m(X) ≥ max(in-degree, out-degree). In-degree is at least the hyperplane-cover number of the
  content; out-degree is at least the same for the reach.
- Define ρ(d) = max over X of content·reach/m. Then R_d ≥ 5824/ρ(d). The lattice is self-dual, so ρ(d) = ρ(9−d).

## 2. Results (all checks PASS)
| level | ρ(d) | registers needed | helpers alive (data = 224) | script |
|---|---|---|---|---|
| 1, 8 | 26 | 224 (+8 centres at dim 8) | 0 / **8** | expSF cert.py |
| 2, 7 | 26 | 224 | — | cert2.py (a), (c) |
| 3, 6 | **24** | **243** | **19 / 19** | e25.py (E25 = 24), rho3.py (quads), rho3t.py (triples ≤ 16 tops) |
| 4, 5 | ≥ 32 (E34 = 32, g(4) = 36) | ≥ 182 | 0 | e34.py, gcap.py |
| 2 (sources) | per-port argument | **336** | **112** | cert2.py (b), (d) |

**New lemma (dim 2, source side), `cert2.py`.** The 52 pairs (s,·) leave ⟨p_s⟩ only through dim-2 out-neighbours.
- Charge m(N) to the ports whose line-edge into N is used: a 2-port node gives m/2 to each of its two ports.
- A charge-1 node reaches at most 26 tops.
- The source's own L2 node reaches at most 24 tops if it is a 2-port node (u ∈ K_s forces the same class), or at most 28 if
  it is a 1-port node. A 1-port node reaching more than 20 tops needs m ≥ 2, because every 6-dim subspace holds ≤ 20 ports.
- So every port needs a charge of at least 3, which gives **R_2 ≥ 336, i.e. ≥ 112 helpers alive at dim 2**.
- The bound is locally attainable for every s, in two ways (`dim2tight.py`):
  - same-class pair node + a 1-port node reaching 32 tops;
  - a 1-port node in K_s + two G26 pair nodes.

**Exact minimum of this per-level relaxation (`helperilp.py`, CP-SAT).** Helpers count as intervals [a,b] alive at dims a..b,
costed f(b−a). The constraints are: ≥ 112 alive at dim 2, ≥ 19 at dim 3, ≥ 19 at dim 6, ≥ 8 at dim 8.
- Minimum **585.3 log-mass**: 112×[2,3], 11×[6,7], 8×[6,8].
- With only the dim-2 and dim-8 constraints the minimum is 514.9.
- This is the exact optimum of the relaxation, and a valid lower bound for the full OR-reach model. It is not the exact
  minimum of the full model.

**Pricing (`price3.py`, same ceiling call as price2.py):**
| helper inventory | ceiling |
|---|---|
| LB inventory: 123 × w1 + 8 × w2 (585) | **×2.031** (was ×2.136 in expSF) |
| dims 2 and 8 only: 120 × w1 | ×2.044 |
| hypothetical hub-tree design: 224 × w7 + 8 × w1 (3,711) | ×1.623 |
| 80 × w7 + 8 × w1 | ×1.919 |

## 3. Answer
- **The stop test is not met.** The lower bound is 585 log-mass, well short of the ≈ 4,500 needed to close the lever at ×1.5.
- It is above the ≈ 500 "go to construction" mark only marginally.
- The lever window narrows to **[×1.36 proved cap for today's inventory, ×2.031 proved best possible]**.
- **No construction was found**, so the best known feasible point is still today's word: ≈ 8,250 log-mass, ×1.258 one-shot.
  Task 2 (explicit design, `wordmodel.py` replay) was not reached: the 3 h went into the bounds.

**Structural facts for a construction.** These were found here and are not exploited yet:
- The 14 syndrome classes are exactly the odd vectors of 4-dim spaces X_σ. X_σ^⊥ contains no ports, so class-pure gathers
  reach nothing.
- A port is orthogonal to exactly 4 ports in its own class, 4 in each of 12 other classes, and 0 in the antipodal class
  (`cls2.py`).
- The ports plus the 8 centre normals are 120 vectors: they look like E8 roots mod 2 (O+(8,2) non-singular points), with
  port-orthogonality equal to odd inner product (52 + 4 = 56). This is suggested by the counts and was not checked.

**First obstruction to a cheap design.** The tight levels are dims 1, 2, 7, 8 (ρ = 26) and 3, 6 (ρ = 24).
- At dim 2 every port needs a 1-port node (sibling not a port) or a third register.
- In a binary "hub-tree" design, each port sits in two 4-dim hubs, and hub contents are crossed with dual hubs at dim 4→5 by
  rectangles X ⊥ W (≤ 36 pairs each).
- A 6-port hub offers only 2 non-port siblings, so this forces ≥ 56 hubs, i.e. 448 registers, i.e. ≈ 224 width-7 helpers
  (≈ 3,700, ×1.62).
- That design is a count-feasible blueprint only. The rectangle cover and the source/target assignment are not verified.

## 4. Next experiment (ranked)
1. **Hub-tree construction.**
   - Gather side: choose 4-dim hubs, each with a flag e ⊂ ⟨e,e'⟩ ⊂ E_X that fixes its pair/quad tree.
   - Source sibling e ∈ K_s; one non-port sibling per port.
   - Scatter side: the dual hubs. Middle: a CP-SAT cover of the 5,824 ordered orthogonal pairs by rectangles ports(X) × ports(W)
     with X ⊥ W.
   - Assign registers to source/target/helper roles, verify OR-reach, then replay literally with `wordmodel.py`.
   - **Stop test:** if no cover exists with ≤ 270 width-7-equivalent helpers (≈ 4,500 log-mass), the lever is very likely
     below ×1.5 in practice. If a cover exists, the relaxed lever is ≥ ×1.5, and the next obstruction is lifting it
     (time order, cancellation, the centre frame growth that the relaxation ignores).
2. Strengthen the lower bound with a joint dims 1–3 LP. Use node types (content, reach, in, out) and flow conservation between
   levels, which this scout's per-level bound ignores. Only a jump of about ×8 (to ≈ 4,500) would close the lever, so this
   is low priority.

New files: cert2.py, e25.py, e34.py, gcap.py, rho3.py, rho3t.py, dim1.py, dim2tight.py, cls2.py, helperilp.py(.json), price3.py.
Note: one `pgrep -f` kill matched this scout's own shell (exit 144). It did not hit any other scout's process.
