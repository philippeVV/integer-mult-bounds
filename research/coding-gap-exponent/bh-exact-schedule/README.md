### First exact helper-free schedule of the Braverman–He code (q = 9 Hermitian Tanner code on PG(2,9))

This is a schedule of BH's 166-session code on the full PG(2,9), with **no helper ports and no dropped sessions**:
- **Gates.** There are 182 gates, each a [10,5] MDS code over F₉. Points use doubly-extended RS rows. Lines use the √q-twisted rows.
- **Order.** The order is Singer cyclic: points P_i at 2(16i mod 91), and lines L_l at 2((16l − 3) mod 91) + 1.
- **Orientation.** It is Eulerian, 5 in and 5 out at every gate. There are |R| = 166 back arcs, which equals dim C, and Γ − R is
  acyclic.
- **Exactness.** A gate-by-gate F₉ simulation finds every gate bijective and **rank(Φ − I) = 0**. So the session map is the identity:
  the schedule is exact and sound, with 166 = S sessions.
  - For comparison, BH's own order needs 9 C ports and 9 J ports.
- **Negative controls.** Rescaling one coordinate of the line codes gives rank(Φ − I) = 1. Rescaling all of them gives 61. These match
  the dense dim W of 165 and 105.

**Reproduce.** You need the file `schedule9.json` (the order, the 910 oriented arcs and R) plus the scripts:

```
python3 check_json.py schedule9.json          # stdlib: perfect difference set, incidences, Eulerian 5/5, Kahn acyclicity, |R| = S(9) = 166
python3 verify.py 3 cyc9.txt dense            # numpy: transportation 166/166, 182 gates bijective, rank(Phi - I) = 0, dense dim W = 166
PERTURB=one python3 verify.py 3 cyc9.txt      # control: rank(Phi - I) = 1
```

`verify.py` rewrites `schedule9.json` in the current directory when the check passes. Our rerun reproduced it byte for byte.

**Scope and caveats.**
- This is a circuit over **F₉**. It is **not** a candidate for `finite_bit_contract.py`.
  - The checker accepts the incidence: 166 roles, 182 gates, and 1076 physical edges (166 + 744 + 166).
  - It rejects the scalar map (`Wrong input-to-output permutation`), because F₉ gates have no XOR form.
  - No binary encoding can supply one: (F₉, +) has exponent 3.
- It is not routable. Girth 6 forces every head→tail path to have ≥ 5 arcs, and 5·166 = 830 > 744.
- **No κ claim.** On this exact, lane-free instance, every frame family we screened is certified at s ≥ Wm, so it gives no saving:
  - GRS-F₂, 23 phase-frame families, m = 4..8, 3 draws: 345 of 345 ranks full;
  - the two-class theorem, 41 partitions: all m ≤ 15;
  - the rational ray test, all scalings: 184 rays, every gcd trivial.
- The schedule is offered as a test bed. It does not open a route.

Prepared with Claude Code (Claude Opus 5.5) assistance; not independently reviewed.
