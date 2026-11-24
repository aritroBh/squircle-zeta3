# HEXC2R — adversarial referee report on HEXC2 (exact level-6 certificate for the hexagon-corner closed form of M₃(2π/3))

**Verdict: PROVED AFTER REPAIRS (computer-assisted).** The mathematical content of HEXC2 is
correct and is verified here end-to-end with fully independent code: the three identities
ID₂, ID₃, ID₄ are exact rational combinations of proved relations among level-6 coloured MZV
values (residual 0 on every word, three independent exact arithmetic backends), through the
Möbius map f(z) = (ζ²z+1)/(ζ⁴z+1). No frame double-shuffle rows are used; quotient dimensions
are 3, 8, 23, 69 against Deligne's 3, 8, 21, 55 (weights 1–4). The repairs are
documentation-level (one sign typo in a displayed identity in HEXC2's result.md), plus two
strengthenings delivered by this referee (exact verification of the pulled-back 1-forms; an
exact bonus certificate for Ψ⁽²⁾(1/2)). The numeric standard of the conversion audit matches
and slightly exceeds the paper's thm:W4R standard.

## The claims, and what I checked

Claims under review: ID₂, ID₃, ID₄ (weights 2,3,4) — hence COR5's closed form
M₃(2π/3) = W₄ + W₃ + W₂ — are exact consequences of proved relations among level-6 CMZVs
(regularised double shuffle, Zhao distribution for d ∈ {2,3,6}, weight-1 branch, lifts)
through the Möbius pullback; residual 0 on every word; no frame double-shuffle rows;
quotients 3, 8, 23, 69 vs Deligne 3, 8, 21, 55.

**Every claim is confirmed**, with the whole pipeline re-implemented independently
(`reflib.py`, ~19k chars, no imports from HEXC2/HEXC1 code).

### Task 1 — the certificate algebra: CONFIRMED (exact + numeric)
- My own shuffle/stuffle/regularisation/distribution generators reproduce cmzv6's relation
  rows **exactly as sets** at weights 1–4 (4/57/682/7869 primitive rows; family counts match).
  (One cache-mutation bug was found *in my own* re-implementation during development and fixed;
  cmzv6's `reg` correctly copies before popping.) [step_01]
- With my own exact code I rebuilt rows and targets and re-ran the certificate:
  Σ yᵢ·row(prowᵢ) − T = **0 on every word** at weights 2, 3, 4 (cert rows 28/229/1695, target
  supports 26/177/920) — with FLINT fmpq_mat (all weights), with an independent pure-Fraction
  backend (weights 2,3), and with my own modular arithmetic mod 5 primes (weight 4). The
  Euler bonus cert (ζ(2) = π²/6 in the weight-2 span) also passes exactly. Max y numerator:
  143 digits, as stated. [step_02]
- Row validity: every row used by the certificates (**all 1952**) plus 141 unused rows,
  classified by first-emit family (RDS, distribution d=2/3/6, weight-1 branch, lifts),
  evaluates to 0 at 65 dps with **my own CMZV evaluator** (worst 1.37e-63). The evaluator is
  written from the definition as truncated iterated sums (geometric tails from a 1/2-split,
  no special functions), cross-validated against both HEXC2 evaluators and against known
  constants (worst 9.21e-63 on 136 words). [step_03, step_04]
- Certificate family usage: w2: 20 DS + 5 dist2 + 1 dist3 + 2 lift1; w3: 203 DS + 11 dist2
  + 2 dist3 + 13 lift1; w4: 1612 DS + 26 dist2 + 5 dist3 + 50 lift1 + 2 lift2. **No frame
  rows.** [step_05]

### Task 2 — the conversion: CONFIRMED (exactly equal to mine; numeric to 70 digits)
- My own conversion (Möbius pullback + Chen path composition through 0 + shuffle
  regularisation + cross-ratio frames), derived and coded independently, produces rational
  vectors **exactly identical** to hexc2_common's on all 38 HPL words and all 9 Li atoms
  (with the same frames). [step_06]
- Numeric audit at 70 dps with my evaluator: all 38 HPL conversions worst 6.19e-71; all 9
  Li atoms worst 1.15e-70; 8 elementary atoms worst 3.33e-64. With the **wrong**
  composition order 36/38 words fail (worst 2.29) — the audit genuinely discriminates the
  order convention (the bug HEXC2 fixed). [step_06]
- The 12 Möbius maps re-enumerated exactly over ℚ(ζ) (my own cyclotomic arithmetic): exactly
  12, exponent tuples identical to mobius_maps.json; task-map images (0,1,2,5,4) exact; each
  root omitted exactly twice. [step_05]
- Cross-ratios of {0, 1/2, 1, −1, ∞} (my own exact rational enumeration): exactly the 15
  values {−3,−2,−1,−1/2,−1/3,1/4,1/3,1/2,2/3,3/4,4/3,3/2,2,3,4}; the closed-form arguments
  {1/3, −1/3, 2/3, 1/2, −1/2, −2} all occur (1/3: 6 frames, −2: 6 frames, ...); **1/9 does
  not**, and Li₄(1/9) = 8(Li₄(1/3)+Li₄(−1/3)) holds (exact even/odd series split; 1.12e-201
  at 200 dps). [step_05]
- Frame-independence of the conclusion: with MY OWN frame choices the targets differ from
  HEXC2's as vectors but vanish numerically (1.78e-66 / 2.92e-65 / 1.82e-64 at 65 dps), and
  the differences lie **exactly** in the proved-relation span (weights 2,3: exact rank
  argument; weight 4: mod 3 primes). [step_07]

### Task 3 — what is proved, what is numeric
**Exact computation:** the certificate algebra (Σ y·rows = T on every word; three backends);
my row generators = cmzv6's (set equality); my conversion = their conversion (rational
vector equality); the 36 pulled-back 1-form identities (exact over ℚ(ζ), this referee's
step_10); the 12 maps and cross-ratio tables (exact); quotient ranks mod 2 primes; the Ψ⁽²⁾
bonus cert (step_11). **Validated only numerically:** (i) that the exact rational vectors Tₙ
built by the conversion code represent the analytic identities IDₙ (the conversion *code* —
validated to ≥ 60–70 digits on every conversion by now three independent evaluators); (ii)
the inputs: the G1 formula M₃ = −ζ(3)Φ + ζ(2)Ψ + R₃ (claimed proved in the paper's notes;
outside HEXC2's scope) and the COR5 closed-form coefficients (lattice-reduction outputs,
validated to 1115 digits here).
**Same standard as thm:W4R?** Yes — the paper's proof of thm:W4R states the conversion code
was "audited and checked numerically ... but did not verify formally". HEXC2's audit is on
*every* conversion (not a sample), with two independent evaluators (three with this
referee's); and the analytic input of the conversion is now verified exactly (step_10), so
the standard is matched and slightly exceeded.
**Could the conversion be made exact cheaply?** Essentially yes: (1) the pulled-back forms
are exact rational 1-form identities (verified, step_10); (2) path composition, rescaling,
reversal and shuffle regularisation are Chen/IKZ theorems — the only un-formalised part is
~100 lines of exact-rational bookkeeping code; (3) a further cheap exact certificate is the
vanishing of the symbol of Tₙ (necessary, exact linear algebra); (4) no Borel/denominator
argument is needed anywhere because π⁴ drops out exactly.

### Task 4 — the numbers: CONFIRMED (1094 digits)
With my own polylogs (direct series for |z| < 1, worst 1.63e-1200 vs mpmath; Jonquière
inversion for Li₄(−2), |diff| 1.24e-1200 vs mpmath) I recomputed the closed form at
1200 dps: M₃(2π/3) = 0.0107173629322894654184095762912248958557189757384643357601556…
|closed form − COR5's 1095-digit value| = 1.03e-1095 (1094-digit agreement, the reference's
own precision). Components W₂/W₃/W₄/F₄/Ψ⁽²⁾ match COR5's hex_*.txt to ≤ 4.8e-1117. The
fully explicit form (with certified Ψ⁽²⁾) also agrees to 1094 digits. [step_08, step_11]

### Quotient dimensions: CONFIRMED
My own vectorised GF(p) elimination (two primes): ranks 2/28/229/1695 at weights 1–4, hence
quotients **3, 8, 23, 69** vs Deligne's conjectured 3, 8, 21, 55. Weights 1–2 are complete;
2 relations missing at weight 3, 14 at weight 4 — the certificate uses none of them.
[step_09]

## Repairs

- **R1 (required, documentation):** HEXC2's result.md displays ID₂ as
  "−H₀,₁(½) − H₀,₋₁(½) + ζ(2)/2 = …". The certified (and numerically correct) identity has
  **+H₀,₋₁(½)** (the fit's coefficient): −H₀,₁(½) + H₀,₋₁(½) + ζ(2)/2 = ½log²2 + Li₂(1/3) +
  ½log²(3/2). The printed form is off by 0.897. The certificate and cert_w23.json are
  correct; only the display in result.md (and any paper sentence quoting it) must be fixed.
- **R2 (clarification):** the exact statements ID₂/ID₃/ID₄ take as inputs (i) the G1 HPL
  formula for R₃, Ψ (proved in the paper's notes, not re-proved here), and (ii) COR5's
  closed-form coefficients (lattice-reduction outputs, here validated to 1115 digits and now
  certified through the same machinery).
- **R3 (optional strengthening, delivered):** the pulled-back 1-form identities are now
  verified exactly over ℚ(ζ) (step_10); the full classical form additionally needs
  Ψ⁽²⁾(1/2), for which I add an exact certificate (step_11, cert_psi2.json):
  Ψ⁽²⁾(1/2) = (3/32)π² − (9/8)Li₂(1/3) + (3/4)log 2 log 3 − (9/16)log²3
  (28 weight-2 rows, residual 0 on all 36 words, multiplier heights ≤ 4 digits).

## The sentence the paper can carry (mirroring thm:W4R)

> **Theorem (hexagon corner, computer-assisted).** M₃(2π/3) = W₄ + W₃ + W₂ with
> W₄ = F₄ + ζ(2)Ψ⁽²⁾(1/2) + (3/8)ζ(3)log 3, where
> F₄ = −21Li₄(−2) − (65/8)Li₄(1/3) + (13/2)Li₄(2/3) + (137/6)Li₄(1/2) + (39/64)Li₄(1/9)
> plus products of lower-weight periods over ℤ[1/6], W₃ and W₂ are in lower-weight periods,
> and Ψ⁽²⁾(1/2) = (3/32)π² − (9/8)Li₂(1/3) + (3/4)log 2 log 3 − (9/16)log²3; in particular
> M₃ at the hexagon corner is classical. Each weight-homogeneous identity is an exact
> rational combination of proved relations among level-6 coloured MZV values (regularised
> double shuffle, Zhao's distribution relations of orders 2, 3, 6, the weight-one branch
> identities, and shuffle lifts), through the Möbius map f(z) = (ζ²z+1)/(ζ⁴z+1) sending
> {0, 1/2, 1, −1, ∞} onto sixth roots of unity; the residual is exactly 0 on every word
> (36/252/1764 words of weights 2/3/4), the multipliers come from exact rational solves
> (heights ≤ 143 digits), π⁴ comes out exactly, and no relations beyond the standard ones
> enter (the 2 weight-three and 14 weight-four relations missing from the standard set play
> no role). The conversion code is audited numerically to ≥ 60 digits on every conversion
> with independent evaluators; its analytic input — the pulled-back differential forms — is
> verified exactly over ℚ(ζ).

## Evidence table (all numbers printed by my code)

| check | step | label | residual / result |
|---|---|---|---|
| my relation generators == cmzv6 row sets, weights 1–4 | 01 | exact | 4/57/682/7869 rows, sets equal |
| cert Σ y·rows = T on every word, w2/w3/w4 + Euler | 02 | **exact** | residual 0 on 36/252/1764/36 words; fmpq + pure-Fraction + mod 5 primes |
| my CMZV evaluator vs numeric_cmzv / fast_cmzv | 03 | cross-check | worst 5.01e-63 / 9.21e-63 (136 words, 65 dps) |
| all 1952 cert rows + 141 unused vanish numerically | 04 | numeric (60+ dps) | worst 1.37e-63, all families |
| 12 Möbius maps, task-map images, omitted roots | 05 | exact over ℚ(ζ) | identical to mobius_maps.json |
| cross-ratio values; {1/3,−1/3,2/3,1/2,−1/2,−2} present; 1/9 absent | 05 | exact rationals | 15 values |
| Li₄(1/9) = 8(Li₄(1/3)+Li₄(−1/3)) | 05 | exact series + numeric | 1.12e-201 @200 dps |
| weight-1 branch rows vanish | 05 | numeric | 3.36e-66, 2.37e-66 |
| my conversion == hexc2 (38 HPL + 9 Li) | 06 | **exact** | rational vectors identical |
| HPL/Li conversions numeric audit | 06 | numeric 70 dps | worst 6.19e-71 / 1.15e-70 |
| elementary atoms (log2, log3, log32, iπ, ζ2, ζ3, ±π², π⁴) | 06 | numeric 65 dps | worst 3.33e-64 |
| wrong composition order fails | 06 | numeric | 36/38 words off (worst 2.29) |
| my-frame targets vanish; differences in proved span | 07 | numeric + exact | 1.78e-66/2.92e-65/1.82e-64; rank equality (exact w2/3, mod p w4) |
| ID₂ display sign typo in HEXC2 result.md | 07 | numeric | printed form off by 0.897 |
| M₃(2π/3) from closed form (own Li₂/Li₃/Li₄ + Jonquière) | 08 | numeric 1200 dps | 1.03e-1095 vs COR5 (1094 digits); components ≤ 4.8e-1117 |
| quotient dims 3, 8, 23, 69 vs Deligne 3, 8, 21, 55 | 09 | exact mod 2 primes | ranks 2/28/229/1695 |
| 36 pulled-back 1-form identities | 10 | **exact over ℚ(ζ)** | all hold |
| Ψ⁽²⁾(1/2) closed form certificate | 11 | **exact** | residual 0 on 36 words, heights ≤ 4 |

## Files (new_results/HEXC2R_referee_certificate)

- `reflib.py` — the referee's independent library (exact word algebra + relation generators,
  truncated-iterated-sum CMZV/HPL evaluator, conversion machinery).
- `steps/step_01.py` … `steps/step_11.py` — the checks above (each appends to progress.log).
- `verify_hexc2r.py` — the < 5 min aggregate verifier (PART A exact; PARTS B–E numeric/exact).
- `cert_psi2.json` — the referee's bonus exact certificate for Ψ⁽²⁾(1/2).
- `result.md`, `progress.log`.
