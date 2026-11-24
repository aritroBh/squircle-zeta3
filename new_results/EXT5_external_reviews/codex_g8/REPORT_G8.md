# REPORT G8: M_2, M_3, M_4 at every angle (Theorem 20, Proposition 5, Theorem 21 of the full version)

Files (all in `paper/ieee/proofs/`): `sec_g8_M2alpha.tex`, `sec_g8_M3alpha.tex`, `sec_g8_M4alpha.tex`, `macros_g8.tex`, `wrap_g8.tex`,
`scripts/g8_*.py` with `.out`, this report. Nothing else was edited. Scripts are run from `scripts/` (they read their data files by relative path).

## Build

`tectonic -C --keep-logs wrap_g8.tex` (own wrapper; it inputs `sec_g6_flat`, `sec_g6_field`, `sec_g6_ladder`, `sec_g6_fourdilog` first, then the three new sections).
- 45 pages in all (21 for the introduction and the four g6 sections, 24 for g8 and the bibliography). New sections: M2 about 10 pages (p. 22 to 32), M3 about 7.5 (p. 32 to 39), M4 about 5.5 (p. 40 to 44).
- `wrap_g8.log`: 0 undefined references, 0 overfull boxes (none at all), 6 underfull boxes (2 in `sec_g6_field`, 4 in the new sections, all small); `pdftotext`: 0 occurrences of "??". PDF (17:36:56) newer than every source (latest 17:36:43).
- Build trap, avoided: the g6 ladder section uses `\mathsf`; after it, math in a `\section`/`\subsection` heading or at a new font size (`\small`, `\footnotesize`) makes tectonic `-C` look for `cmss12` at that size, which is not in the cache (error "Font OT1/cmss/m/n/14.4 not loadable"). The new sections therefore have no math in headings and no size changes. A master that prints g8 after the g6 ladder will hit the same trap if math is added to those headings.

## Hand / machine / not proved, per result

### Theorem 20 (M_2 at every angle)
- Hand: reduction to the rational part (Lemma g6:prop:formula, Lemma g6:lem:first(b)); the criterion that an identity in the algebra of nested sums gives an identity for every n (Lemmas g8:lem:formal, shift, hcoef, criterion); `W_3 = M_2 - W_2` is the combination `Y(c)+Y(-c) - zeta(2)(1-c^2)/(2c) artanh c`; method, ranges, branches and domains for the classical form.
- Machine: M1+M2 `g8_allN.py 2` (0.3 s, re-run yes): all-n identity for the rational part, odd coefficients, coefficients c^N with N <= 1 (word-by-word equality in Q(n)). M3 `g8_cor3_verify.py` (18 s, re-run yes; copy of the notes' COR3k script; 33/33 items): 6 weight-two HPL identities, 5 dilog reductions, duplication, 10 weight-three identities by exact differentiation (derivative identity, value at 0+, ranges). M4 (same script, "TASK 2"): substitution; all 31 coefficient keys in Q(c) agree with the classical form; zeta(3) coefficient 7/4 + 2c.
- Not proved: the Euler-substitution/IBP integral representation of W_3 (notes G1, Theorem C); uniqueness; any relation among 1, pi^2, zeta(3). Two gap boxes in the section.
- The classical form is stated as a Corollary (`g8:thm:classical`), not a Theorem, to avoid a clash with the paper's numbering.

### Proposition 5 (M_3 at other angles)
- Hand: definitions (V, symbol of H_w, B_2 = Lambda^2 V, delta_{2,2}); delta_{2,2} kills shuffle products and Li_4 (Lemma g8:lem:vanish); one term cannot: exhaustive over every single term r(c)G with r in Q(c)^x and G of any rational-coefficient symbol in any letters (constants projected out): the two components of delta_{2,2}(Omega) have ratio (1-c)/(1+c), not in Q (Lemma g8:prop:one); the four differentials of Li_{a,b} (Lemma g8:lem:diff); specialisation at x = 1/2; square-corner collapse.
- Machine: M5 `g8_allN.py 3` (0.6 s): M_3 = R_3 - zeta(3)Phi + zeta(2)Psi for all n (defines Omega through the 74-term R_3; the notes had `b3_allN_proof.py`). M6 `g8_M3_delta22.py` (<0.1 s, re-run yes): S(Omega) = 24 pure tensors; delta_{2,2}(Omega) = (1-1/c)(x^y)^(y^z) - (1+1/c)(x^z)^(y^z), equal to the paper's formula; controls (900 random and all 243 shuffle products, 16 Li_4 symbols: 0; positive control nonzero). M7 (same script): generic delta_{2,2}(Li_{2,2}) from the 15-tensor two-variable symbol, specialised: delta_{2,2}(two terms) = delta_{2,2}(Omega) exactly, constant letters removed.
- Established, precisely: (a) the formula and non-vanishing for the explicit 56-term expression; (b) impossibility for one term over the stated class; (c) the two-term identity at symbol level modulo the constant letter [2].
- Not proved: invariance of delta_{2,2}(Omega) under the choice of expression, hence "M_3 is not classical as a function" (needs independence of HPLs and classical functions over Q(c)); Goncharov's converse; anything at one angle (period conjecture, hexagon, c_0 = 1/3, golden corners, corner lists); the notes' "no combination of depth-two MPLs at S-unit arguments repairs the [2]-residual" (COR8B; not reconfirmed by COR8BR); that two terms plus classical functions equal Omega as functions. One gap box, 7 items. The impossibility "with constants" that the referee reports flagged is not claimed.

### Theorem 21 (M_4 at every angle)
- Hand: the zeta(3) part through Lemma g6:thm:ladderprecise (z = zeta(3), gamma_1 = 1): [s3]P_{4,n} = -[s2]P_{3,n}; the pi^4 part from two rational functions by partial fractions and three elementary sums; the all-n argument written as mathematics (Lemma g8:prop:crit and four points in the section).
- Machine: M8 `g8_allN.py 4` (3.5 s, re-run yes): parts 1 and zeta(2) (109 and 22 terms; 29 and 5 words; even identity, odd identity, c^N for N <= 1). M9 (same run and `g8_pi4.py`): [s4]P_{4,n} = (n+3)/((2n-1)(2n+1)(2n+3)), [s2^2]P_{4,n} = -n/((2n-1)(2n+1)(2n+3)); ladder rung as an identity in the algebra.
- Not proved: the consequences the paper lists after the theorem (depth-two statements at generic angle, hexagon, c_0 = 1/3, golden corners); uniqueness of the 109 + 22 terms; the real-number split of the parts (the theorem is proved with formal coefficients). One gap box, 4 items.

## New argument, needs referee (file:line)

1. All-n proof of Theorem 20 for k = 2 (the paper says "proved for every n by reducing both sides to nested harmonic sums"; the notes have it exact only through c^400, G1): `sec_g8_M2alpha.tex:113-121` (M1/M2 box and output), criterion `:41-93`. Same method and library as the k = 3, 4 cases refereed by G2R, but nobody has read this application.
2. The mu side computed from Definition g6:def:formal for symbolic n inside `g8_allN.py` (the notes used a formula `mu4_formula.pkl` checked for n <= 640 and then converted): `sec_g8_M4alpha.tex:57`, `sec_g8_M2alpha.tex:113`. Tested against `g6_common.py` for n = 1..8 (all monomials, k = 2, 3, 4) and by mutants.
3. Shifts S_j with j >= -2 and the dictionary at N = 2n + delta, delta >= -1 (`sec_g8_M2alpha.tex:49-75`); the validity range is used to show no denominator vanishes for n >= 1 (`:93`). The notes' referee audited the denominators, not this formulation.
4. Lemma g8:prop:one (one term cannot, exhaustive class) `sec_g8_M3alpha.tex:142-156`; Lemma g8:lem:vanish `:129-140`.
5. Hand derivation of the four differentials of Li_{a,b}, `sec_g8_M3alpha.tex:158-186`, and the specialisation to the constant-free iterated-differential symbol, `:194-209` (the notes verified the differentials numerically; COR8BR by finite differences).
6. Ranges/branch table for the eight arguments, `sec_g8_M2alpha.tex:130-135` (the notes had a grid check).
7. Square-corner collapse as a computation, `sec_g8_M3alpha.tex:125`.
8. The zeta(3) and pi^4 parts of Theorem 21 by hand (the notes had an all-n identity for the explicit 15-term Psi): `sec_g8_M4alpha.tex:55,61-71`.

## Disagreements between `main.tex` (draft 47) and the notes

1. Thm 20 proof sentence claims an all-n proof; the notes (G1) have an exact check to c^400 and say an all-n proof was not carried out. Now supplied (item 1 above). The paper's "120 digits" for the classical form rests on COR3k values; we re-checked 40 digits.
2. Thm 20 statement "W_3 is an iterated integral of weight three whose letters are all real": proved here only as a Q(c)-combination of weight-three HPLs and zeta(2) times a weight-one HPL, letters 0, 1, -1. The integral representation (G1 Thm C, Euler substitution) is not reproved.
3. Thm 21 "22 such terms of weights one to three": the file has 21 HPL terms plus one constant term 1/2 (3, 6, 12 terms of weight 1, 2, 3). R_4 has 1/8/27/73 terms of weight 2/3/4/5 (the notes' G2 result.md says 3/8/28/70, a documentation error found by G2R; total 109 correct). R_3 has 2 + 16 + 56 = 74 terms; Psi has 14 terms and a constant.
4. Prop 5, two terms: delta_{2,2}(Omega) = delta_{2,2}(two Li_{2,2}) holds in the constant-free symbol (modulo the letter [2]). With [2] kept as a formal letter, four components remain (COR8B); our script prints them. The text says "(of the function, modulo constants)" for delta_{2,2}(Omega) but not for the two-term sentence; add "modulo constants" there. "Checked by two independent programs" is right (COR8B, COR8BR, plus ours); COR8R did not confirm it (convention: sign reversed, shuffle-inconsistent constants; report "partly verified").
5. COR8B's printed closed forms for delta_{2,2}(Li_{2,2}), delta_{2,2}(Li_{3,1}) had an error (found by COR8BR; harmless at x = 1/2); our generic delta_{2,2}(Li_{2,2}) agrees with the corrected formula.
6. The paper's numbers (Theorem 20, Proposition 5, Theorem 21) are used for the three statements only; every auxiliary result is a lemma or corollary numbered inside its section.

## Gap boxes (all)

- M2: (i) Euler-substitution representation of W_3 not reproved; (ii) real-number identity only, no uniqueness or independence claims; (iii) how Y was found plays no role. Second box: machine steps are not formally verified programs; inputs (eleven terms of Y, the table of Lemma wt3) are part of what is checked.
- M3: one box, seven items (Goncharov's criterion quoted; invariance under the expression; numbers at one angle and corners; one term only for constant-coefficient G; [2]-kept residual and the notes' impossibility not claimed; differentials proved by hand; the 74 terms are data).
- M4: one box, four items (consequences not treated; terms are data; formal vs real split; numerics are not proofs).

## Checks (fresh, own scripts)

- `g8_numeric.py` (7 s): M_2 and M_4 from the definition (Cauchy integral, 40 digits) against Theorem 20 (HPL form), the classical form, and Theorem 21 at alpha = pi/3, pi/2, 2pi/3, 1.8: differences 2.3e-39 to 4.5e-39 (M_2), 1.4e-37 to 1.8e-37 (M_4); M_4(pi/2) = 0.00602799311260912250... as in the paper.
- `g8_allN.py selftest`: 98 words x 7 values of n x 6 shifts = 4116 checks of the HPL-coefficient rewriting against the definition; harmonic numbers; stuffle. `g8_mutants.py`: 4 mutants of Y and 2 of the k = 4 data are all caught.
- mu side vs `g6_common.py` at n = 1..8 for every monomial (k = 2, 3, 4). `g8_pi4.py`: partial fractions, closed forms, 1/2880 combination (exact).
- Reports: G1 verify.md (flat formulas), COR3k (planner re-run of the ten identities: 34 items; ours 33), G2R (all-n method for k = 4, mutation tested, verified), COR4/COR8R/COR8B/COR8BR for Prop 5 as described in the section. They confirmed the notes, not this text.

## Bibliography

No new entry needed. Cited: Goncharov95, GGSVV, Gangl16 (in `refs.bib`), Lewin81, Zagier, DLMF, Gon98, Hof00, RV, Rud76, sympy, mpmath, flint, arb (all in `refs.bib` or `supplement.tex`); `wrap_g8.tex` lists them.
