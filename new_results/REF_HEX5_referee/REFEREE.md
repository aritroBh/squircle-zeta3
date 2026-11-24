# Referee report: prop:order5, thm:hexagon5, prop:octagon5depth (main.tex draft 34, lines ~1686-1738)

All code here is the referee's own. No project code is imported. The project's data are only read:
`K1g_hexagon_order5_assembly/coef_acc_raw.pkl` (the 307-monomial table) and
`OCT5D_octagon_order5_depth/parsed_mons.pkl` (the octagon integrand). The project's verifier
`OCT5D2_octagon_depth_test/verify_oct5d2.py` was re-run once, read-only (no bytecode, no new files).
Python: `.venv/bin/python` (mpmath, sympy).

## Verdicts

- **prop:order5: PROVED AFTER REPAIRS (wording only).** The law is re-derived exactly by independent code.
  Every printed number checks. Repairs R1, R2, R6.
- **thm:hexagon5: PROVED AFTER REPAIRS (computer-assisted).** No coefficient is fitted. The table matches
  the referee's direct computation to about 110 digits on six hexagons. Repairs R1 to R5. One statement is
  wrong as printed: the denominator is L^9, not L^5 (R3).
- **prop:octagon5depth: PROVED at the symbol level** (computer-assisted). The pairings were reproduced
  exactly by an independent computation (own integrand, own symbol code). The number-level statement is
  conditional on the period conjecture, as the paper says. Wording repairs R7 to R9.

## 1. Statement checks (numerical, independent)

**Direct method (`direct5.py`).** Polar perimeter integral
P(p) = ∫ F^{-1/p} sqrt(1+(F_θ/(pF))²) dθ, F = Σ_j w_j^{-d} l_j^{1+d}, l_j = c_j|cross(V_j,u)|, d = p-1.
This is exactly the weighted family w_j|k_j⟨v_j,x⟩|^p with w_j k_j = c_j.
- Differentiation under the integral sign: at each node the integrand is a truncated power series in d
  (orders 0..5). The k-th coefficient times k! is integrated.
- Own tanh-sinh, all orders at once, step halving until two levels agree. Each sector is split at its
  midpoint and each half is parametrised from its own corner, so the vanishing term is |V_j| sin φ, with no
  cancellation.
- Finite differences in p (`fd5.py`, no differentiation under the integral) check the interchange.

**Closed form (`formula5.py`, `hpl.py`).** The 307 coefficients were read from `coef_acc_raw.pkl` and
lambdified. E_e = 120·Σ coef·(HPL monomial)·x^k·ζ's, x = 1-r. Edge data (A0, A1, B, r, q0, L², |V0|²,
weights) are computed from scratch from the corners. The 13 HPLs are the referee's own: series at 0, then
analytic continuation along the real axis by Taylor re-expansion of the ODE system. Self-test: Li_n,
-log(1-x), S_{1,2} and Li_{3,2}(x,1) match mpmath or nested sums to 1e-50 at x = -7.5..0.8
(`log_hpl_selftest.txt`). On r = 1 edges the monomials have poles, so E_e is taken as the average of
x = ±10^{-k} at raised precision (error O(10^{-2k})).

| hexagon (corners V0,V1,V2; weights) | ratios r | P^(5)(1) (direct = formula) | digits (dps 110) |
|---|---|---|---|
| regular, w = 1,1,1 | 1, 1, 1 | 656.41536196762475447664760905831120964826426… | 109.8 |
| paper example (1,0),(2,2),(0,1), w = 1,3,1 | 2, 1/2, 1 | -449.23396800163978089849025140306005176022… | 110.7 |
| generic A: (1,1),(-27/10,33/10),(-2,1), w = 5,1,2 | 13/10, 1/2, 20/13 | -1474.6784800498412080730263501189100249635… | 110.2 |
| generic B: (1,0),(4,16/5),(0,1), w = 1,2,3 | 4, 5/16, 4/5 | -1448.7675291780296112161919413483521378596… | 109.9 |
| golden (1,0),(φ,1),(0,1), w = 2,3,5 | φ, 1, 1/φ | -78.358599873964602127531124178168617182820… | 109.5 |
| generic C: (2,-1/3),(26/7,3/14),(1,1), w = 3/2,4,1 | 3/2, 7/5, 10/21 | 720.74638522850829292847185026162593439825… | 109.9 |

- Each edge also agrees separately, to the same precision (`log_hex5_dps60.txt`, `log_hex5_dps110.txt`).
  Residuals are at working precision (1e-58 at dps 60, 1e-108 at dps 110).
- The r = 4 edge (x = -3) and the r = 1 edges (limit x → 0) agree as well as the others.
- Regular hexagon: the direct value equals the paper's 656.4153619676247544766476… and HEX5's 113 digits.
  Golden hexagon: it equals the paper's -78.3585998739646… and the notes' 55-digit formula value.
- Finite differences (`log_fd5_generic_B_r4.txt`, hexagon B, central 5th difference, h = 10^{-3}/2^k,
  six Richardson levels): -1448.76752917802961121619194135, which is 5.0e-30 from the jet value. So the
  interchange of d/dp and the integral is confirmed.
- Side observation: steps h ≥ 0.025 do not converge (`log_fd5_generic_B_r4_largeh.txt`). This fits the
  paper's own result that the diamond-end series has radius zero. P is C^∞ but not analytic at p = 1, so
  P^(5)(1) must be justified by differentiating under the integral, not by a series argument (R4).

**prop:order5 on fully local polygons (`local_examples.py`, `log_local_examples.txt`).**
- An exact P^(5) per edge (own exact derivation, `exact_local_edge.py`) equals the direct value to 1e-58
  on every edge. Its ζ(5) and ζ(2)ζ(3) coefficients equal the printed Z_e^(5) and Z_e^(2,3) exactly.
  - Parallelogram ±(3,1), ±(1,2): Z^(5) = 0 on both edges, because the foot is at a corner and h = L.
    It is a weak test, so the next two cases were added.
  - Paper's sheared parallelogram ±(2,0), ±(1,2): ζ5 = -6912√5/25 - 172800√13/28561 and
    ζ2ζ3 = -10752√5/25 - 821760√13/28561 (as in the notes).
  - Affine-regular, non-regular hexagon (1,0),(1,1),(0,1), weights 1,2,4: ζ5 = 135 + 675√2/32 and
    ζ2ζ3 = 15 + 195√2/32.
- Integer relations from the direct numbers alone (no exact code) give the same relation at two precisions:
  - parallelogram (3,1),(1,2): 6 constants, dps 110 and 150;
  - affine-regular hexagon, weights 1,2,4: 17 constants including log 2, dps 160 and 200, maxcoeff up to
    4.8e5, residuals 1e-156 and 1e-197.
  The ζ5 and ζ2ζ3 coordinates equal Z_e/L exactly on every edge.
- Regular 2m-gons, m = 2..10: Σ_e Z_e from the printed law equals the printed C_m^(5) and C_m^(2,3) to
  1e-47. The values for m = 3 are 2835√3/32 and 1395√3/32; for m = 2 they are -90√2 and -210√2.
- Diamond: the exact P^(5)(1) = √2(-90ζ5 - 210ζ2ζ3 + 396ζ2² + 720ζ3 - 480ζ2), as printed
  (`log_misc_checks.txt`).
- Regular hexagon: the exact P^(5)(1) minus the printed 15-term closed form is **0** (sympy). P'''(1)
  matches as well, and P'(1) = 2√3 log 2 (`log_exact_regular_hexagon.txt`).

**r = 2 example (`pslq_r2.py`, `log_pslq_r2.txt`).** The r = 2 edge of the paper hexagon was computed in
the weighted convention, not the notes' unit normals, from the referee's own jets at dps 250 and 320.
Integer relation on the 20 constants:
- Both precisions give the same relation, maxcoeff 921600, residuals 3.5e-246 and 8.2e-316.
- [ζ5] = -20187/2000 √5 and [ζ2ζ3] = -858/125 √5, as printed.
- Li4(1/2), Li5(1/2) and Li4(1/2)log 2 occur, with coefficients -1242/25, 396/125 and 594/25 (times √5).
- The lower-weight coordinates differ from E5e's unit-normal relation (for example 96 log 2 against
  -153 log 2). This confirms the paper's remark that the convention changes other terms but not ζ5 or
  ζ2ζ3.
- The local law on this edge gives -2079/250 √5 and -477/250 √5 (`log_misc_checks.txt`), so the
  non-local parts -711/400 √5 and -1239/250 √5 check.
- The 13 HPLs at x = -1 (own code) agree with the exact reductions quoted in K1g/RUL1 to 3.6e-61.

## 2. Derivation check

| step | how established in the notes | referee |
|---|---|---|
| prop:order5: only (a,b) with a+b = 5, a,b ≥ 1 carry ζ5; only (3,2), (2,3) carry ζ2ζ3; [ζ5] = 24B, [ζ2ζ3] = -12B | exact (E5 step_03), numeric table 4.8e-35 | Re-derived by hand (Faà di Bruno on log B, ψ^(m) at integers), and exact sympy over 84 cases, 0 mismatches (`beta_rule.py`) |
| prop:order5: Z_e^(5), Z_e^(2,3) | exact computer algebra (E5 step_04/05) | **Re-derived independently** (`local_law.py`). At order d^5 a log-degree-5 monomial is leading-log, so its coefficient is [s0^a s1^b] of \|d/dt(y/F)\| with F = 1 + A0 t(e^{s0}-1) + A1(1-t)(e^{s1}-1). It depends only on A0 = ξ0, A1 = ξ1 and the geometry, so it does not depend on the weights or the far lines. Derived minus printed = 0 symbolically, for both laws and for the parallelogram short forms. |
| C_m^(5), C_m^(2,3) | exact, symbolic in h | numerical, m = 2..10, 1e-47 |
| regular hexagon closed form | found first by PSLQ (HEX5, 125/145 digits); then exact (Codex; RUL1 target 3) | **Re-derived exactly** by independent code (`exact_local_edge.py`): the full d^5 integrand of a κ = 0 edge, integrated with Beta derivatives. Difference from the printed form = 0. |
| thm:hexagon5: E5 integrand | exact; its d³ part equals K1d's E3 (80 monomials) | not re-derived symbolically; confirmed through the final values |
| J-table (c ≥ 1) | exact engine: series in x, e-sums, and an exact shuffle reduction to the 13 Lyndon HPLs; quadrature check 1.3e-76 | not re-derived symbolically. The atom reductions are exact algebra (step_25: shuffle algebra over Lyndon words; no PSLQ). |
| J-table for r ≥ 2 | **not argued**: the x-series converges only for \|x\| < 1, that is 0 < r < 2 | Needs R2. Confirmed numerically at r = 4 (x = -3) and on the r = 2 edge. |
| assembly into 307 monomials | exact sympy (K1g step_30/34) | Not re-run symbolically. Checked numerically on 18 edges, 110 digits. |
| standalone ζ5, ζ2ζ3 coefficients = Z_e for every r | exact at 6 random points (K1g step_58) | Re-done against the referee's derived law at 4 more random rational points: 0 (`zeta_coeff_check.py`). Only the x^0 monomials carry a standalone ζ5 or ζ2ζ3. |
| poles x^{-6}..x^{-1} cancel | "structural", numerical κ → 0 | Immediate: each J(n;a,b,c)(x) = ∫ t^n log^a t log^b(1-t) log^c(1-xt) dt is analytic for x < 1. Confirmed by the r = 1 limits. |
| denominator "L^5" | K1g says L2^{5/2} | **Wrong.** The table has L2^{-1/2} up to L2^{-9/2}, and 87 coefficients really have an L^{-9} pole (random data, `log_denominators*.txt`). This is dimensionally forced: q0 and \|V0\|² enter to degree 5 at order five. See R3. |
| differentiation under the integral | not stated for order 5 (it is stated for thm:hexagon) | Needs R4. FD confirms (5e-30). |
| r = 2 fold into 20 constants | E5e: PSLQ (336 digits). K1g's first "exact" fold used atom values at x = -1 **found by PSLQ** (K1g step_49/50, `atom_mzv_m1.pkl`). RUL1 later proved all 13 reductions with DS+DUAL certificates. | Certificates not re-checked. Atom values agree numerically (3.6e-61). The relation was re-found independently in another convention. |

**Fitted parts.** No coefficient of thm:hexagon5 or prop:order5 rests on PSLQ. PSLQ was used in three places:
1. The regular-hexagon closed form (HEX5). It was later derived exactly, and the referee re-derived it.
2. The r = 2 relation (E5e).
3. The 13 atom values at x = -1 in K1g's original fold (step_50). These were replaced by RUL1's exact
   certificates, which the paper now cites ("by double shuffle and duality relations (in our notes)").

The folder name `E5e_pslq_340` refers to item 2 only. It is a check, not an ingredient of either result.

**Residual trust.** Three steps were not re-derived symbolically by the referee:
- the weight-5 integrand engine;
- the J-table engine;
- the final assembly.

All three are exact computer algebra in the notes. They are confirmed to about 110 digits on 18 edges with
generic data: r ∈ {5/16, 10/21, 1/2, 1/φ, 4/5, 1, 13/10, 7/5, 3/2, 20/13, φ, 2, 4}, weights from 1
to 5. A wrong coefficient that survived this would have to vanish on all 18 edges.

## 3. prop:octagon5depth

**What is claimed.** The object is the part of total logarithmic degree five of the order-five edge
integrand of P^(5)(1)/u. These are the leading-log terms, including log 2 and log(1+√2). Each integral
∫ t^n Π logs is rewritten as multiple polylogarithms G(w;1) with letters 0, 1, 2+√2, -(1+√2). Its symbol
(constant entries included, over the S-units 2, 1+√2) has a non-zero antisymmetrised (3,2) component.
The two pairings are -1485+1035√2 and 900-540√2. Classical polylogarithms and products pair to 0. So, at
the symbol level, no combination of Li5's and products has this symbol. For the number, this needs the
period conjecture, and the paper says so.

**Independent check (all exact, own code).**
- `oct_weight5_fast.py` re-derives the weight-5 integrand from scratch in exact Q(√2) arithmetic. It uses
  own geometry and a leading-log expansion of \|d/dt(y/F)\|. All 252 rows (1303 coefficients) equal
  OCT5D's `parsed_mons.pkl`: **0 mismatches**.
- `own_symbol.py` computes the symbol with its own Goncharov formula, I(0; w_n..w_1; 1), dropping zero
  entries, and its own S-unit factorisation and pairing. Two inputs are written by the referee: the
  top-weight formula ∫t^n G(w;t) → (1-w_1^{n+1})/(n+1) G(w;1), and the shuffle of single-letter G's.
- Controls:
  - G(0⁴,p;1) = -Li5(1/p);
  - G(p)G(q) = G(p,q)+G(q,p);
  - Li5 at the 7 S-units pairs to 0 (both entry orders);
  - shuffle products pair to 0.
- Result: local block (0,0); single-p = single-q = (-53955/8+4770√2, -10215/2+14445/4√2);
  mixed (48015/4-8505√2, 11115-15525/2√2); **total (-1485+1035√2, 900-540√2)**, exactly as printed.
  The pairing is the same for either order of the symbol entries.
- The printed correction C = (-3600+2160√2)Li_{3,2}(b,1) + (11880-8280√2)Li_{3,2}(1/√2,1) pairs to the
  same values.
- The project verifier, re-run, prints 21 PASS (`log_rerun_verify_oct5d2.txt`).
- The referee's own jets give the octagon's P^(5)(1) = 164.12775911330597065104567331245436130507711746936…
  This agrees with OCT5D's 120-digit value to 2e-59 (`log_misc_checks.txt`).

**Does the conclusion follow?**
- Yes, it follows that no combination of Li5's and products has this symbol. A non-zero (3,2) cobracket
  excludes classical polylogarithms; that is a theorem.
- That depth two *suffices* is shown at symbol level by C.
- "Needs depth two" for the *number* is conditional on the period conjecture. The paper says so in the
  statement.
- C is *one* correction; it is not unique (OCT5D2 reports a 1-dimensional symbol kernel of Li_{3,2}
  combinations).
- The non-proportionality to M4's (84,-38√2) holds, but it depends on the normalisation of the pairing.

## 4. Hypotheses

- **Non-degenerate hexagon.** As at order three, the hexagon must be strictly convex: 0 < B < 1, r > 0,
  A0 = 1 - Br > 0. A hexagon with a straight angle is a parallelogram (log A0 diverges).
- **Weights.** P^(5)(1) depends on the weights; Z_e^(5) and Z_e^(2,3) do not (they depend only on
  ξ0 = A0 and ξ1 = A1). The theorem refers to the weighted family through "the notation of
  Theorem hexagon". Weights from 1 to 5, including non-integers, were tested.
- **Range of r.** The J-table is a power series in x = 1-r, valid for |x| < 1. r ≥ 2 needs analytic
  continuation (R2). The tests cover r = 2 and r = 4.
- **r = 1.** The monomials have poles at x = 0, so at r = 1 the formula is a limit (R3).
- **Differentiation under the integral.** For |p-1| < 1/4 the p-derivatives up to order five are
  dominated by C(1+|log φ|^5)φ^{-1/4} near each corner, which is integrable. FD confirms (5e-30). The
  paper does not say it at order five (R4).
- **Independence.** "Coefficient of ζ(5)" is unique only modulo the conjectural linear independence of
  the constants. This applies to prop:order5 ("all of the ζ(5)" of a fully local polygon) and to the
  r = 2 numbers, where the paper does say it.

## 5. Repairs and sentences that claim more than is established

- **R1 (prop:order5, statement).** "The local terms of an edge $e$ contribute exactly
  $Z_e^{(5)}\zeta(5)+Z_e^{(2,3)}\zeta(2)\zeta(3)$ to $P^{(5)}(1)$" and "If the polygon is fully local,
  these are all of the $\zeta(5)$ and $\zeta(2)\zeta(3)$ of $P^{(5)}(1)$."
  → "Evaluated with derivatives of the Beta function, the local terms of an edge $e$ give
  $Z_e^{(5)}\zeta(5)+Z_e^{(2,3)}\zeta(2)\zeta(3)$ plus products of $\zeta(2)$, $\zeta(3)$, rationals and
  logarithms of the corner term values. If the polygon is fully local, this is all of the $\zeta(5)$ and
  $\zeta(2)\zeta(3)$ of $P^{(5)}(1)$ in this presentation (as coefficients of the number, assuming the
  constants are linearly independent)."
- **R2 (thm:hexagon5, proof).** Add: "The table of integrals is derived from power series in $x$, which
  converge for $|x|<1$. Both sides are analytic in $x<1$, so the formula holds for every $r>0$; we
  checked $r=4$."
- **R3 (thm:hexagon5, statement).** "with coefficients that are polynomials in $A_1$, $B$, $r$, $q_0$,
  $|V_0|^2$, $L^2$ and $1/x$, divided by $L^5$ (the poles at $x=0$ cancel)."
  → "… divided by $L^9$ (the poles at $x=0$ cancel; at $r=1$, $E^{(5)}_e$ is the limit as $x\to0$)."
  Also "For every centrally symmetric convex hexagon" → "For every strictly convex centrally symmetric
  hexagon and all weights $w_j>0$".
- **R4 (thm:hexagon5 and prop:order5, proofs).** Add one sentence: "We may differentiate under the
  integral sign, because for $|p-1|<\frac14$ the first five derivatives in $p$ are dominated by an
  integrable function." The same sentence appears for Theorem hexagon. It matters because P is not
  analytic at p = 1.
- **R5 (notes / reproducibility).** The paper points to "our notes" for the 307 coefficients. They exist
  only as a 4 MB sympy pickle (`coef_acc_raw.pkl`), and the notes say "denominator L2^{5/2}", which is
  wrong. Add a plain-text table to the supplement and correct the denominator in K1g/result.md.
  The table needs E_e = 120·Σ, x = -κ = 1-r, lA_j = log(A_j/w_j), A0 kept as a symbol equal to 1-Br.
- **R6 (prop:order5, proof paragraph, optional).** "Exact symbolic integration shows that the law minus
  the raw integrals is exactly 0" is correct. An independent derivation (leading-log expansion) now gives
  the same law; it can be cited as a second derivation.
- **R7 (prop:octagon5depth, statement).** "The weight-five part of $P^{(5)}(1)/u$ has a nonzero $(3,2)$
  part."
  → "Write $P^{(5)}(1)/u$ as the integral of its order-five edge integrand, and let $W$ be the terms of
  total logarithmic degree five, written with multiple polylogarithms $G(w;1)$ with letters
  $0,1,2+\sqrt2,-(1+\sqrt2)$. The symbol of $W$ has a nonzero antisymmetrised $(3,2)$ part."
- **R8 (paragraph after the proposition).** "At the level of symbols the weight-five part of
  $P^{(5)}(1)/u$ is classical plus the exact correction $C$ … but the remaining constants are not yet
  identified."
  → "… is classical plus an exact correction $C$ (one choice; it is not unique)". Delete "but the
  remaining constants are not yet identified": the last paragraph of the section says the closed form is
  now known (numerically).
- **R9.** "depth two thus enters with the octagon's non-local terms at order five, as in the squircle at
  $M_4$" → add "(at the level of symbols)".
- **Table II.** "Thm. hexagon, Thm. hexagon5 … 3 checks each" can add "referee: six hexagons,
  110 digits". The row "Prop. octagon5depth … Proved† (symbol level)" is accurate.

## Files

- `direct5.py`: polar jets, any centrally symmetric polygon, orders 0..5; scalar P(p)
- `hpl.py`: own HPLs (`log_hpl_selftest.txt`)
- `formula5.py`, `run_hex5.py`: the table vs direct (`log_hex5_dps60.txt`, `log_hex5_dps110.txt`)
- `fd5.py`: finite differences (`log_fd5_generic_B_r4.txt`, `log_fd5_generic_B_r4_largeh.txt`)
- `beta_rule.py`: Beta rule (`log_beta_rule.txt`)
- `local_law.py`: independent derivation of Z_e^(5), Z_e^(2,3) (`log_local_law.txt`, `local_law_derived.pkl`)
- `exact_local_edge.py`: exact κ = 0 edges; regular hexagon (`log_exact_regular_hexagon.txt`, `exact_regular_hexagon.pkl`)
- `local_examples.py`: fully local polygons, PSLQ, regular 2m-gon sums (`log_local_examples.txt`)
- `zeta_coeff_check.py`: standalone ζ5, ζ2ζ3 of the table (`log_zeta_coeff_check.txt`)
- `denominators_check.py`, `denominators_check2.py`: the L^9 denominators (`log_denominators.txt`, `log_denominators2.txt`)
- `pslq_r2.py`: r = 2 relation (`log_pslq_r2.txt`)
- `misc_checks.py`: octagon value, r = 2 local law, diamond, HPLs at -1 (`log_misc_checks.txt`)
- `oct_weight5_fast.py`: octagon weight-5 integrand vs OCT5D (`log_oct_weight5_fast.txt`, `own_weight5_rows.pkl`)
- `own_symbol.py`: own symbols and pairings (`log_own_symbol.txt`)
- `log_rerun_verify_oct5d2.txt`: re-run of the project verifier
- `t_regular.py`: quick sanity run of `direct5.py`
