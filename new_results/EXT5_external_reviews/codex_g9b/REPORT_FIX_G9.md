# REPORT_FIX_G9: fix pass on Proposition 9, Theorem 26, Proposition 10

Backup of the three sections, macros and wrapper: `backup_before_fix_g9/`. Build `wrap_g9.tex`: 58 pages, 0 Overfull boxes, 0 undefined references or citations, no "??" in the PDF text, no em dash in the three sections. Table `scripts/g9_hex5_table.txt` unchanged, sha256 f0fd1f7a771301273e14f4f8f94a5e109e50f78c10c3265a23dc080bf9c2030a (verified); the 15-minute chain was not re-run.

## Statements
The three quotes follow the approved draft-49 wording (Prop 9 clause "plus products of zeta(2) and zeta(3) with polynomials in the edge data and the logarithms of the corner term values"; Thm 26 final parenthesis; Prop 10 "Its symbol has a nonzero antisymmetrised (3,2) part: in the pairing ..." and "at any arguments"). A diff against `paper/ieee/main.tex` (draft 47) shows only these places and cross-reference text. The definitions of xibar, xi_+, xi_x were inside the quoted Prop 9 statement; they are now in the setting paragraph, so the quote is verbatim. A "What we prove" remark follows each quote.

## Findings
1. Lemma z5, factor 120 (hex5, proof of Lemma z5): bracket coefficients are 24 int sum g_ab and -12 int (g32+g23); times 120 they give the local law. Corrected. New script `g9_hex5_z5factor.py` (exact; also reads the two table lines): 120 x coefficient - law = 0 for both, pass, 1.5 s.
2. W5 defined twice (oct5): T_5 is the whole integral, W_5 the weight-five part, defined by the explicit formula with A(w,n). All later uses consistent.
3. Lemma topwt: restated with the limit lim (tau^{n+1}-1)G(w;tau) = 0 for w_1 = 1; A(w,n) = 0 then. Remainder has weight at most k-1.
4. Split-entry rule stated; cited `g9_oct5_symbol.py` lines 104-121 (113, 115-116, 118-119). New script `g9_oct5_coincident.py`: Sym G(a,a;1) = 1/2 Sym G(a;1)^2 for a = p, q, and G(a)G(b) = G(a,b)+G(b,a) for five pairs; pass, 0.03 s.
5. Gap (i) removed. New Lemma `g9:lem:arb` (projection pi: U -> V, pi|V = id, U = C^x tensor Q): Pi = pi^{tensor 5} fixes Sym W5, commutes with shuffles, sends Li5 symbols to the shape killed by Lemma zero (a). Lemma zero extended to coefficients in C, to both orders of the entries, and to arbitrary tensors in (b). Independence of [2], [1+sqrt2] proved by the norm. The proof of "no combination of Li5 at any arguments and products" is now complete by hand, no citation needed. I checked the argument myself, including the reversed-order case.
6. W not closed under the product: now the graded algebra W = Q + sum W_l. `g9_hex5_lyndon_det.py` (executes the matrix construction of `g9_hex5_s2_atoms.py`, 0.3 s) prints determinants 1, -2, -6, -96, 5760 (the referee's values).
7. Thm 26 last sentence: "What we prove" paragraph: statement about the formal representative; for the number it assumes linear independence of the constants (as Prop 9). The approved paper wording now carries it too.
8. Golden hexagon: item 7 of the text-after list rewritten. Not reproduced: 3e-32. Reproduced: finite differences 3.9e-29 (`g9_direct5.py`, `g9_direct5_golden.out`); referee jets, exact derivatives, table limit vs jets 1.1e-51 (`referee_g9/r9_direct.py`, `r9_misc.py`, `r9_misc.out`).
9. `scripts/g9_planner_check_beta.py` copied unchanged, run: diamond equal to M2; sheared z5 = (-197413632 sqrt5 - 4320000 sqrt13)/714025 = -6912 sqrt5/25 - 172800 sqrt13/28561, the other four coefficients agree; 1.9 s, `g9_planner_check_beta.out`. Cited in Checks of Prop 9; the "no script" remark is withdrawn (gap note (iii) rewritten). The writer's disagreement 4 in `REPORT_G9.md` (read only) is withdrawn here.
10. r = 2: symmetry (x,y) -> (y,x) maps the r = 2 edge to the r = 1/2 edge (isometry preserving K_p and the weights 1,3,1); both give -83.5391363475376327634345072141820934358 (`referee_g9/r9_r_half.out`). Stated plainly that the r = 2 corollary rests on the notes' RUL1 certificates, not re-derived (also in the Not-proved list).
11. Unit-normal parenthesis: two-line argument added (item 5 of the text-after list): the z5 and z2z3 lines have no l_j; a shift l_j -> l_j + c log 2 only produces log^i 2 (i >= 1) times weight <= 4 constants, never the monomials zeta(5), zeta(2)zeta(3).
12. Minors: 1302 non-zero octagon coefficients everywhere (1303 in the script output is the union of keys with a notes table that stores one explicit zero; explained in P1; verified 1302 from the rows); oct5 em dash removed; pole list 0, -1, ..., -(n+1), 1/(N+2) dropped; harmonic range i <= 5 for the local table (H^(5) occurs at a+b = 5, e.g. (0,5)), i <= 4 for the non-local sums; "genuine multiple polylogarithms" replaced (H00101 = Li_{3,2}(x,1), H01011 = Li_{2,2,1}(x,1,1) not of the classical form; no minimality claim); Codex 7 (S2 label in the "by hand" list; continuation phrased on the punctured slit plane with removable singularity at 0); Codex 8 (occurrence is a statement about the presentation); Codex 10 (double shuffle justified by finite-cutoff stuffle and limit, since the word (-1) is conditionally convergent); Kimi F3 (chain time machine dependent: 896/899/901 s). Kimi F2 void (g4_edge.py exists). Opus F9 (exact J check): only cited as the referee's result, the script needs the chain's jtab5.pkl, not kept. Opus F16: referee jets cited as the sharper check. "Which report confirmed" paragraphs updated in all three sections.

## Scripts re-run
g9_hex5_z5factor.py pass 1.5 s; g9_hex5_lyndon_det.py pass 0.3 s; g9_planner_check_beta.py pass 1.9 s; g9_oct5_coincident.py pass 0.03 s; g9_oct5_rows.py (0 mismatches, 0.9 s); g9_oct5_symbol.py (totals -1485+1035 sqrt2, 900-540 sqrt2, 5.6 s). New outputs next to the scripts.

## For the planner (suggested paper changes beyond the three approved)
- Golden hexagon text after Thm 26: replace "finite differences ... to 3e-32" by "an independent computation of P^(5)(1) from the polar perimeter integral (derivatives in p taken exactly under the integral sign) to 1e-45, and finite differences to 4e-29".
- "H00101, H01011 are not [Nielsen], so at order five genuine multiple polylogarithms appear" -> "are not of this form, so harmonic polylogarithms of depth two and three occur in the table".
- Paragraph after Prop 10: "no polylogarithm (Theorem 26)" -> "(Proposition 9 and the closed form above; Theorem 26 at r = 1)".
- Table II row for Prop 10: label "Octagon: nonzero (3,2) component of the weight-five symbol" (Codex 9).
- Optionally add to the text after Thm 26 that the table is also checked exactly against its defining series through x^30 (referee).
- "= 120 sqrt2 a_5" (Prop 9) is not proved in the supplement (gap note (ii)); keep or drop.
- "independent finite differences agree to 2e-8" (octagon paragraph) is far weaker than what is now available (3.6e-31 here).

## Unsure
- W_5 is defined by the procedure (shuffle expansion plus Lemma topwt), not intrinsically; the text says so.
- The split-entry rule is taken as the definition and tested by script; no proof that it equals the symbol of the tangentially regularised integral is written out.
- The new text of all three sections has not been read by an independent pass.
