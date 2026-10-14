# Leaf-destination transport: the exact maximum is max(1, floor(log2 n))

**This is a scoped written proof, checked by its author with exact computation. Nobody has reviewed it independently yet. It makes no κ, precision-invoice or rounding claim.**

Card: `uniform-transfer-audit` (track adversarial-audit, priority revision 4, chooser nonce `9d450e657d814dfba70c4f6537e46fac`; see `TASK-CARD.md`). It addresses the open quantified question recorded on main in two places:
- `research/leaf-destination-transport`, which asks whether `max(1, floor(log2 n))` holds for all n;
- PR1416's `single-register-interval-transport`, which asks to "determine the exact maximum between h−1 and h+2".

## Result

Take PR628's balanced dirty word with free, noncommuting, contractive factors, in both orientations. Inject a unit error at the destination of a completed leaf, and transport it to any later cut. Then:
- **every register has free-algebra mass at most h = ceil(log2 n) for n ≥ 2**, and E(1) = 1;
- together with PR1416's family achieving h − 1, **the exact maximum is h − 1 or h**. The previous range was h − 1 to h + 2.

The proof is in `PROOF.md`. It rests on three lemmas:
- **Lemma S:** the destination row's source coefficient is 0 or s·P_K.
- **Lemma L2:** every entry of every prefix map is 0 or a single signed monomial.
- **Lemma L1:** each entry of the exit column has mass at most 1.

A 12-case table over the lowest common call then finishes the argument. Its single case with both a dirty source and a dirty destination cancels exactly by Lemma S.

## Follow-up: the bound is exact, E(n) = max(1, floor(log2 n)) for every n (`PROOF2.md`)

- **Upper bound for all n.** Two lemmas give it:
  - Lemma X (exit support): after an injection inside a call of size k, at most floor(log2 k)+1 registers are nonzero when the call ends, and at most floor(log2 k) of those are scratch.
  - Lemma I: a sign − call is the exact leaf-by-leaf inverse of the sign + call.
  - In every case where the source is dirty, the exporting phase is the floor half, so the count is floor(log2 m).
- **Attainment for every n ≥ 2.** An explicit family achieves the bound. With top sign +, inject at the f_0 leaf along the chain U, U′, U, U′, …, and read out at the first f_{n−1} leaf. The exact value there is
  `y = Σ_{l=1}^{H} (−1)^{H−l} f_{n−1}…f_{n−⌈n/2^l⌉} · f_{⌊n/2^l⌋−1}…f_1`,
  whose distinct tail lengths prevent any cancellation, so the mass is H = floor(log2 n).
- **Exact checks** (`check2.py MODE LO HI [famHI]`):
  - n = 1..80, both orientations and all injection/cut pairs: Lemma I, L1, Lemma X (the bound and its tightness), the per-case closed forms, E(n), Lemma P and the family. See `check_ok_*.log`.
  - Family up to n = 512; closed form for n = 2..512 (`fam_closed.py`).
  - Controls: no reversal gives E(7) = 54; a ceil-half split gives E(6) = 3; dropping the restoration breaks Lemma I (`controls2.log`).
- **Scope.**
  - The lower bound is proved for the (U,V) case with sign +. The other source-dirty cases match the census but are not proved.
  - The proof has not been reviewed independently.

## Exact checks

| Script | What it checks |
|---|---|
| `census.py N...` | Exact census of the envelope. The values equal floor(log2 n) for n = 1..128; n = 33..128 is new, and the logs for n = 71..128 are `ext_*.log`. |
| `chainform.py` | A second, independent method. It agrees with `census.py`. |
| `PYTHONPATH=. python3 lemmas.py 64` | Lemmas L1, L2 and S hold exactly for n ≤ 64, both signs. Runs in about 53 s. |
| `controls.py` | Mutation controls. Not reversing the negative phases, or dropping the restoring phase, breaks the lemmas; the envelopes become 3/55/32,691,955 and 2/11/291 instead of 2/3/4. |
| `lcabound.py`, `patterns.py` | With independent entries in place of Lemma S's paired entries, the envelope becomes h + 1. So Lemma S is necessary. |

The bound h is attained at n = 2^h for h ≤ 7 and is strict otherwise.

## Open next questions

- Cover multiple or accumulated errors, internal leaf values, rounding and aliasing.

Credit: woklerx (leaf-destination-transport), cbethin (PR1416), Arham Tahir (PR1374), Hanson Sun (PR1085), Chafik Boukhalfa (PR628). Prepared by Philippe Van Velzen with substantial Claude (Anthropic) assistance. A sub-agent wrote the proof and the census; the main session re-ran the census for n ≤ 64, the lemma checks and the controls. Apache-2.0.
