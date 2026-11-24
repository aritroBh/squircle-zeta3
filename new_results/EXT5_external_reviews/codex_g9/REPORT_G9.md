# REPORT G9: order five on polygons (Proposition 9, Theorem 26, Proposition 10 of the full version)

Writer: group G9. Date: 2026-11-23. Statements copied from `paper/ieee/main.tex` (draft 47; numbers from `main.aux`). None of the three results is in the core paper: `paper/ieee/core/main.tex` only says "\fullref{prop:order5}, \fullref{thm:hexagon5} and \fullref{prop:octagon5depth} of the full version" (l. 1268, l. 1410). No existing file was edited; the `*g8*` files were not touched.

## Build

`cd paper/ieee/proofs && tectonic -C --keep-logs wrap_g9.tex` (final build, after the last edit of any `.tex` file):

| item | value |
|---|---|
| pages | 53 (pp. 1-27: the five g4 polygon sections, unchanged; pp. 28-51: the three new sections; pp. 52-53: references) |
| Proposition 9, `sec_g9_order5.tex` | pp. 28-34 (7 pages) |
| Theorem 26, `sec_g9_hex5.tex` | pp. 35-45 (11 pages) |
| Proposition 10, `sec_g9_oct5.tex` | pp. 46-51 (6 pages) |
| undefined references / citations | 0 (`grep -c undefined wrap_g9.log`) |
| overfull boxes | 0 (`grep -c Overfull wrap_g9.log`; none above 2 pt either) |
| `??` in `pdftotext` | 0 |
| PDF newer than sources | yes |

Warnings left: `Underfull \hbox`, and hyperref "Token not allowed in a PDF string" from section titles of group g4 (as in `wrap_g4.log`).

## Files

New in `paper/ieee/proofs/`: `sec_g9_order5.tex` (Prop 9), `sec_g9_hex5.tex` (Thm 26), `sec_g9_oct5.tex` (Prop 10), `macros_g9.tex` (only `\providecommand`), `wrap_g9.tex` (preamble of `wrap_g7.tex`; inputs `sec_g4_poly`, `sec_g4_edgelaw`, `sec_g4_Cm`, `sec_g4_hexagon`, `sec_g4_nonlocal` first so that every reference resolves), `wrap_g9.pdf`, `wrap_g9.log`, this report. New in `scripts/`: the `g9_*` files (scripts and their `.out`). The table of the 307 coefficients is `scripts/g9_hex5_table.txt` (5 348 394 bytes, sha256 `f0fd1f7a771301273e14f4f8f94a5e109e50f78c10c3265a23dc080bf9c2030a`). Labels start with `g9:`. The large intermediate pickles of the chain (4 MB each) were written to a scratch directory and are not kept; `G9WORK=<empty dir> sh scripts/g9_hex5_chain.sh` re-creates them (896 s).

## Results

### Proposition 9 (order five, local terms of an edge), pp. 28-34

Inventory category B. After this section it is **by hand, except for two finite exact machine identities**.

*By hand.* Lemma `g9:lem:beta5` (weight-five Beta moments: weight bound; `[z5] = 24 B` iff a+b=5, a,b>=1; `[z2 z3] = -12 B` iff (a,b) in {(3,2),(2,3)}); Lemma `g9:lem:loc5` (which monomials carry z5, z2z3; `Z5 = 2880 int sum g_ab`, `Z23 = -1440 int (g32+g23)`); Lemma `g9:lem:psi5`(a,b) (the series Psi on a local edge as one explicit power series in two variables, the weight-five version of `g4:lem:mixed`); the formal meaning of "contribute ... plus products of zeta(2), zeta(3), rationals and logarithms"; fully local polygons; the parallelogram short forms (substitute xi0 = xi1 = 1); the regular 2m-gon sums `C_m^(5)`, `C_m^(2,3)` (hand computation from the edge data of `g4:thm:Cm`).
*Machine steps (exact, re-run yes).* M1 `g9_edgelaw5.py`: the polynomial identity of Lemma `g9:lem:psi5`(c) in Q[xi0,xi1,u0,m] (series to degree 5, 21 monomials, integrate in t): both differences exactly 0 (0.5 s); the same script reproduces Z_e of Prop 7 at order three (0) and the 2m-gon forms with symbolic s (0). M2 `g9_local5.py`: exact P^(5)(1) of the diamond, the regular hexagon and the sheared parallelogram +-(2,0), +-(1,2) from the full length element (4 s); all equal the paper's / notes' values exactly (difference 0); a regular-hexagon edge has 91 monomials in the d^5 integrand (paper: 91). M3 `g9_beta5.py`: 504 formal Beta moments, rules checked, quadrature 4.8e-35 (9 s).
*Claim R7 (regular hexagon closed form):* **exact, not PSLQ.** M2 integrates the whole d^5 integrand of the fully local hexagon exactly and gives the 15-term form, with 2835 sqrt3/32 and 1395 sqrt3/32. PSLQ (HEX5, 125 and 145 digits) only found it first; the referee's `exact_local_edge.py` derived it independently.
*Numerical only.* `g9_direct5.py` (own code; perimeter of K_p from the definition, finite differences with Richardson, 100 digits, h0 = 0.002, 7 levels): regular hexagon 1.7e-30, sheared parallelogram 6.3e-28 against the closed forms.
*Not proved (gap box).* Independence of the constants (statement about the number); the equality with a_5 for the diamond (the paper's Proposition 1); the notes' "second exact derivation" and the 100-digit quadrature (not repeated); the machine steps are not written out by hand.

### Theorem 26 (P^(5)(1) of every hexagon), pp. 35-45

Inventory category B. **Is the 307-monomial closed form derived exactly end to end? Yes, by exact computer algebra (rational arithmetic in SymPy), with no integer-relation search in any of the 307 coefficients.** Where PSLQ appears in the notes:
1. `HEX5` regular hexagon: PSLQ found it; it is derived exactly (Prop 9, M2). Not an ingredient of Thm 26.
2. `E5e_pslq_340` (relation at r = 2, 250/280/310/340 digits, maximal coefficient 2 574 432 000, LLL cross-check): a **numerical identification**, not an ingredient. The paper's r = 2 coefficients (-20187/2000 sqrt5 and -858/125 sqrt5) are reproduced **exactly** by folding the 307-table at x = -1 with the exact values of the 13 functions (`g9_hex5_r2fold.py`: all 20 coefficients of the referee's weighted-convention relation agree; edge value -83.5391363475376327634345072141820934358). The folder name refers to this item only.
3. K1g's first values of the 13 functions at x = -1 were found by PSLQ; RUL1 replaced them by exact certificates (finite double shuffle and duality, rational arithmetic). We did not re-derive the certificates; we re-ran the notes' verifier in place (41/41 PASS, 57 s, `g9_rul1_verify.out`) and checked the 13 values numerically at 110 digits (9.2e-111). This concerns only the r = 2 corollary.

What is not independent: the J-table engine and the assembly (S4-S6) are one implementation (the notes' code), re-run by us (896 s) with exact equality of all 307 coefficients with the notes' table (0 differing). The weight-five integrand engine (S1) is confirmed by a second code (Theorem 25's `g4_edge.py`): symbolically (`g9_hex5_enginecheck.py`, 336 monomials at order 5 and 0 differing coefficients at every order 0..5, 967 s, compared with the notes' `E5_polys.pkl`, read only) and at three rational points of the edge data (630 coefficients each, 0 differences, 42 s). The label "Proved† (computer-assisted; refereed)" in Table II is therefore accurate; I neither upgrade nor downgrade it. As for Theorem 25, the weakest links are the J-table and the assembly: exact computer algebra confirmed numerically (J-table 1.3e-76; closed form against the definition 1e-29 here, 1e-110 by the referee), not rederived by hand or by a second code.

*By hand.* First sentence (reduction, differentiation under the integral sign: `g4:lem:smooth` holds for every k); shape of E_5 and the denominator **L^9** (Lemma `g9:lem:shape5`; the notes' body says L2^{5/2}, main.tex draft 47 already says L^9); the series for J (Lemma `g9:lem:Jseries5`; analytic in x outside [1,infinity), J = O(x^c)); facts on HPLs (Lemma `g9:lem:hpl`: MPL = HPL, shuffle, stuffle, and the 13 Lyndon atoms generate every H_w, |w| <= 5, by an exact invertible 16x16 matrix, with no use of linear independence of functions); the continuation to every r > 0 and the cancellation of the poles at x = 0 (Lemma `g9:lem:precise5`); Lemmas `g9:lem:nocst`, `g9:lem:z5` (the coefficients of zeta(5), zeta(2)zeta(3) are exactly Z_e^(5), Z_e^(2,3) for every r; new hand proof, uses Prop 9).
*Machine steps (exact, re-run yes).* S1 E_5 engine: 336 monomials (126 local + 210 non-local), 170 s. S2 26 harmonic sums to the 13 Lyndon generators, 0 s. S3 local Beta table symbolic in N, 40 s. S4 J-table, 210 closed forms, 555 s (+160 quadrature checks, worst 1.3e-76). S5-S6 assembly into 307 monomials, 66 s. Compare, census and table, 45 s. Chain total 896 s. Census: x-exponents -6..0 (39, 44, 44, 44, 44, 44, 48 monomials); all 13 generators occur; H00101 and H01011 in 30 J-keys each; z-exponents (a,b,c): (0,0,0) 210, (1,0,0) 55, (0,1,0) 27, (2,0,0) 13, (0,0,1) 1, (1,1,0) 1; powers L^-9..L^1, odd; largest weight 5. M6: 120 [z5] and 120 [z2 z3] minus the laws of Prop 9 = exactly 0.
*Numerical (fresh, own code).* `g9_hex5_numcheck.py` (the table with series-summed HPLs, 60 digits) against `g9_direct5.py` (definition + finite differences): hexagon A (r = 1.3, 0.5, 20/13): 8.2e-31; hexagon C (r = 1.5, 1.4, 10/21): 6.2e-29; against the referee's 110-digit values: 7e-39 and 7e-40 (limit of our 60 digits). `g9_hex5_r1limit.py`: table at x = +-1e-12 (150 digits) against the exact local value of an r = 1 edge: 2.0e-22. `g9_direct5.py` on the golden hexagon (-78.358599873964602127531..., 3.9e-29 from the notes) and the paper's hexagon (-449.2339680016397809..., 3.4e-31). Notes' atom check re-run: 78 checks, 2.9e-79; e-sum engine: 15 checks, 1.1e-79.
*Not proved (gap boxes, 3).* (1) S1/S4/S5/S6 are exact computer algebra, not written out line by line, J-table not independently re-implemented; (2) linear independence of the 307 functions; (3) the 13 values at x = -1 rest on RUL1's certificates. Not re-run: the notes' 336-digit relation, the unit-normal conversion of the r = 2 edge (1e-55 in the paper), the referee's six-hexagon comparison (two hexagons redone with our own code). Hexagons with some r >= 2 are not covered by our numerical check (the referee covered r = 4); for them the continuation argument is by hand (Lemma `g9:lem:precise5`).

### Proposition 10 (depth two at order five, regular octagon), pp. 46-51

Inventory category B. **A statement about symbols only**; the section says so in its first paragraph and in the gap box.

*Fixed in the text:* the edge and its four terms; the symbol of G(w;1) by Goncharov's recursion with entries +-2^alpha (1+sqrt2)^beta written alpha X + beta Y, symbols in Q(sqrt2) tensor V^(x)5; the weight-five part W_5 (leading-log part of E_5, shuffle expansion, top weight of int t^n G); the (3,2) pairing with Lie elements (copied from Section IX-E); "non-zero" means in Q(sqrt2), exact arithmetic with fractions.
*By hand.* Lemma `g9:lem:topwt` (top weight of int t^n G); Lemma `g9:lem:zero` (the pairing kills Li_5(x) when x, 1-x are S-units, and shuffle products; finite hand checks for the three Lie elements); non-proportionality (two elements of Q(sqrt2) differ already in the rational part, -495/28 against 270/19); the logical consequence.
*Machine steps (exact, re-run yes).* P1 `g9_oct5_rows.py` (the referee's code, copied): re-derives the 252 rows (1303 coefficients) of the weight-five integrand from the geometry; 0 mismatches with the notes' table (0.9 s). P2 `g9_oct5_symbol.py` (referee's code): controls (G(0^4,p;1) = -Li5(1/p); the shuffle identity for G(p)G(q); Li5 at seven S-units pairs to 0 in both orders; random shuffle products pair to 0) all true; blocks: local (0,0); single-p = single-q = (-53955/8+4770 sqrt2, -10215/2+14445/4 sqrt2); mixed (48015/4-8505 sqrt2, 11115-15525/2 sqrt2); **total (-1485+1035 sqrt2, 900-540 sqrt2)**, as printed; the paper's correction C pairs to the same pair (5 s). The notes' verifier re-run in place: 21/21 PASS, 2.9 s (`g9_oct5_verify.out`).
*Numerical.* `g9_direct5.py oct`: P^(5)(1) = 164.127759113305970651045673312, 3.6e-31 from the referee's value (the number, not the symbol statement).
*Not proved (gap box).* (i) "No combination of Li5's and products" for Li5 at arbitrary arguments rests on the cited theorem (Matveiakin-Rudenko, `MatRud`; Charlton-Gangl-Radchenko-Rudenko, `CGRR`) that the (3,2) part of the cobracket of classical polylogarithms vanishes, and on our pairing being a component of it; we prove the S-unit case and products only. (ii) The period conjecture. (iii) The symbol identity S(W_5) = S(classical + C) (the notes' 32/32 rank computation), the classical remainder, the numerical identification of the remaining constants, the decagon, the 12-gon, the nine-Li_{3,2} closed form: not part of Prop 10, not reproduced (only: the pairing of C equals the obstruction).

## "New argument, needs referee" (file:line in `paper/ieee/proofs/`)

1. `sec_g9_hex5.tex:97-103` (Lemma `g9:lem:nocst`) and `:126-134` (Lemma `g9:lem:z5`): the non-local integrals contain no z5 and no z2z3, so the coefficients of zeta(5), zeta(2)zeta(3) in the 307-table are exactly Z_e^(5), Z_e^(2,3) for every r. The notes check this exactly at 6 random points (K1g step 58), the referee at 4 more; the proof is ours (machine-confirmed symbolically in `g9_hex5_compare.out`).
2. `sec_g9_hex5.tex:47-56` (Lemma `g9:lem:shape5`): hand proof of the denominator L^9 and of the odd powers of L.
3. `sec_g9_hex5.tex:76-93` (Lemma `g9:lem:hpl`(d)): the 13 Lyndon atoms generate all H_w, |w| <= 5 (shuffle-algebra homomorphism, exact matrix inversion).
4. `sec_g9_hex5.tex:107-124` (Lemma `g9:lem:precise5`): continuation to every r > 0 and pole cancellation for the weight-five table (analogue of `g4:lem:assembly`; referee's repair R2).
5. `sec_g9_order5.tex:44-65` (Lemma `g9:lem:beta5`), `:71-81` (Lemma `g9:lem:loc5`), `:87-108` (Lemma `g9:lem:psi5`): hand reduction of Prop 9 to one explicit series and one polynomial identity (M1). The notes integrated the whole E_5 by computer (E5 steps 01-05).
6. `sec_g9_order5.tex:116-134`: hand computations of the parallelogram short forms and of `C_m^(5)`, `C_m^(2,3)` (the notes: symbolic machine proof in h and numerical checks).
7. `sec_g9_oct5.tex:39-45` (Lemma `g9:lem:topwt`) and `:67-79` (Lemma `g9:lem:zero`): hand proofs of the integration-by-parts formula and of the two vanishing controls (in the notes the controls exist only as scripts). The definition "weight-five part = sum of (1-w_1^{n+1})/(n+1) G(w;1)" is the referee's; the text states it as the meaning of the paper's phrase.
8. The formal meanings of "contribute ... plus products", "in this form" and "weight at most five" (`sec_g9_order5.tex`, "What the words mean"; `sec_g9_hex5.tex`, same title), in the manner of `g4:sec:local`.

## Disagreements between `main.tex` (draft 47) and the notes

1. Thm 26: main.tex says "divided by L^9"; the body of `K1g/result.md` says L2^{5/2} (corrected only in a planner note at its end). main.tex agrees with the referee and with the table.
2. Prop 10: main.tex's statement reads "It has a nonzero (3,2) part", without the referee's repair R7 ("the symbol of W has a nonzero antisymmetrised (3,2) part"). R8 and R9 are in the paragraph after it. The section proves the R7 form and says so.
3. Prop 9: the notes write Z in expanded variables (w, e, q, Q, p, s); main.tex uses xi-bar, xi_+, xi_x, L0L1/L^2, h^2/L^2. Equal: M1 compares main.tex's form with the raw integrals (difference 0), and M6 compares it with the table.
4. Prop 9, proof paragraph: "a second exact derivation (a log-Beta expansion of the arc length)" has no script in the folders I was pointed to (the referee's `exact_local_edge.py` is an independent exact derivation); `g9_local5.py` is a third.
5. Thm 26 text: golden hexagon "agrees with finite differences to 3e-32": our finite differences (100 digits, 7 levels) agree with the notes' value to 3.9e-29 only. Not a disagreement of values.
6. Everything else reproduced: 307 monomials, 13 atoms, 336/210/126, 91 monomials, 76 digits (we get 1.3e-76), 2835 sqrt3/32, 1395 sqrt3/32, -20187/2000 sqrt5, -858/125 sqrt5, -711/400 sqrt5 and -1239/250 sqrt5 (from the exact fold and the local law), 656.4153619676247544766476..., 164.12775911330597065104567331245..., -78.3585998739646..., (-1485+1035 sqrt2, 900-540 sqrt2).

## Gap boxes in the sections

`sec_g9_order5.tex:38` (Prop 9); `sec_g9_hex5.tex:95`, `:136`, `:175` (Thm 26); `sec_g9_oct5.tex:63` (one box, "Not proved in this section"). Computer-assisted boxes: M1, M2 (Prop 9), P1, P2 (Prop 10).

## Bibliography needs

None. In `paper/ieee/refs.bib`: `Deligne10`, `Goncharov01`, `Kolbig86`, `Lewin81`, `MatRud`, `CGRR`, `Zhao10`, `mpmath`, `sympy`. Already cited in `supplement.tex` (not in `refs.bib`): `Chen77`, `Hof00`, `RV`, `IKZ06`, `Racinet02`. `wrap_g9.tex` has its own `thebibliography` with these entries. "Ree's theorem" is named once in `sec_g9_oct5.tex` with a finite hand check for the three Lie elements used, so no entry is needed; the general statement would need one (Reutenauer, *Free Lie Algebras*, or Ree 1958).

## What a referee should rerun

In `scripts/`, with `G9WORK` set to an empty directory: `python g9_beta5.py` (9 s), `g9_edgelaw5.py` (0.5 s), `g9_local5.py` (4 s), `g9_direct5.py {reg6,shear,golden,paper,oct} 0.002 7 100` (3-11 s each), `g9_oct5_rows.py` then `g9_oct5_symbol.py` (1 s, 5 s), `sh g9_hex5_chain.sh` (896 s), then `g9_hex5_compare.py` (45 s), `g9_hex5_numcheck.py {A,C}` (22 s), `g9_hex5_r2fold.py` (22 s), `g9_hex5_enginecheck_pts.py` (42 s; `G9E5=<E5_polys.pkl>`), `g9_hex5_enginecheck.py` (967 s; same variable), `g9_hex5_r1limit.py` (1 min; `G9TABLE=<coef_acc_raw.pkl>`). The notes' verifiers, run in place and read only: `OCT5D2_octagon_depth_test/verify_oct5d2.py` (3 s), `RUL1_level2_exact_reductions/verify_rul1.py` (57 s). The output of every run is in the `.out` file with the same stem.

Note: the scratch directory holding the chain's pickles was emptied by the environment once during the work (not by me). The table and all `.out` files were produced before that, from the regenerated pickles, except `g9_hex5_enginecheck.out`, `g9_hex5_enginecheck_pts.out` and `g9_hex5_r1limit.out`, which read the notes' pickles (`E5_polys.pkl`, `coef_acc_raw.pkl`, read only; the table is equal to the regenerated one by the exact comparison of `g9_hex5_compare.out`).
