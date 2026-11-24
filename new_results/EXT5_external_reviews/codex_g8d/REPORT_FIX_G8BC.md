# REPORT_FIX_G8BC: fix pass for group g8 (M3func, W3int, merged wrapper)

Date 2026-11-23. Backup of every edited file: `backup_before_fix_g8bc/` (sections, macros, the three wrappers, and `scripts/` copies of
`g8b_check.py`, `g8c_structure.py`, `g8c_numeric.py`). New file: `scripts/g8b_pullback.py` (+ `.out`).

## Result
- `wrap_g8.tex` builds all five sections (M2alpha, W3int, M3alpha, M3func, M4alpha) with `macros_g8.tex`, `macros_g8b.tex`, `macros_g8c.tex`.
  `tectonic -C --keep-logs --keep-intermediates wrap_g8.tex`: 68 pages, 0 Overfull boxes, 0 undefined references/citations, no "??" in the PDF text.
- Corollary G8B' is proved in full; no gap forced a retreat to "rational in c". One step of the referee's sketch was wrong as written and is repaired (see "Unsure").

## Findings
**T1 / Opus A1.** `sec_g8_M3func.tex`, new subsection "Arguments rational in the exponential of half the angle" (line 237 on): statement for every non-constant phi in C(t) (angle case phi(z)=(z+1/z)/2; also t=e^{i alpha/(2N)}, t=tan(alpha/4)), six-step proof, machine step P1. Steps: (1) disc B around the accumulation point, C = M o phi as meromorphic functions; (2) H_w o phi in L_{Sigma_t}, length <= |w|, top part s_w phi_*^{(x)|w|} rho(w) (Lemma rec, chain rule); (3) top component of M o phi; (4) delta22 natural in linear maps, pulled-back (2,2) part = (1-1/phi) A - (1+1/phi) B; (5) phi_* injective, so Lambda^2Lambda^2 phi_* injective, coefficients non-zero; (6) contradiction with Lemma g8b:prop:classical and uniqueness.
`scripts/g8b_pullback.py` (P1, exact, 0.8 s, data file g8_terms_k3_part1.json): phi=(z-i)(z+i)/(2z), phi-1=(z-1)^2/(2z), phi+1=(z+1)^2/(2z); [c] -> e_i+e_{-i}-e_0, [1-c] -> 2e_1-e_0, [1+c] -> 2e_{-1}-e_0; rank phi_* = 3; rank Lambda^2Lambda^2 phi_* = 3; pulled-back delta22 computed directly in T^4 from the 24 tensors equals the naturality formula, 21 non-zero components of 45 (as Opus).
Scope edits: theorem title ("not classical for arguments rational in c"), section title, intro paragraph (line 3), the "In words" sentence (old line 38), "What is not claimed" (2) and (3), gap note (line 317), rests-on table (new row). Stated scope: arguments and coefficients rational in e^{i alpha/2} (= rational in cos and sin of alpha/2); not covered: algebraic arguments of positive genus, values at single angles, whether depth two suffices.

**T2 / Opus A2, A5; Codex 7.** Line 23 repaired (no zeros/poles of arguments, no value in [1,inf), principal branches; Li_2(2c-1) example; piecewise expressions excluded via G8B(a) on a subinterval). "One term" paragraph (line 234): h may have poles in the disc. Definition (line 20), Lemma g8b:prop:classical, proof of (a), "One term": "for a chosen finite presentation of C, let E(C) be the union of ..." (E(C) belongs to the presentation).

**T3 / Codex 2.** `sec_g8_W3int.tex` line 32: regularised definition (I(empty)=1, I(0^r)=log^r/r!, recursion otherwise, simplex notation with trailing-zero regularisation, remaining improper integrals converge absolutely). The claim that only (0,0) needs a convention is removed; proof of (b) (line 206) rewritten accordingly.

**T4 / Codex 3.** Statement, proof, checks and gap note now give: 64 integrals of products of two logarithms (10,11,13,15,15 for j=0..4; 23 diagonal, 41 off-diagonal); 105 ordered regularised length-three iterated integrals after shuffle expansion, 21 ending in 0, coefficients in Q+Qc+Q/c; counts of this expansion, no minimality. New step [7] in `scripts/g8c_structure.py` (4.5 s). Output:
```
[7] counts
    products of two logarithms against one form (nonzero coefficient): 64
    by outer index j = 0..4: [10, 11, 13, 15, 15]
    diagonal (a=b): 23 ; off-diagonal (a<b): 41
    distinct ordered words (j,a,b) with nonzero coefficient: 105   = diagonal + 2*off-diagonal = 105
    uncollected summands (diagonal counted twice): 128
    words ending in 0 (innermost letter 0): 21
    words (j,0,0): [(1, 0, 0), (2, 0, 0), (3, 0, 0), (4, 0, 0)]
    the pure word (0,0,0) occurs: False
    every coefficient times c is a polynomial in c of degree <= 2 (so the coefficients lie in Q + Qc + Q/c): True
```
My counts equal Codex's on every number.

**T5 / Opus B1 = Codex 6; B3.** W3int line 180: "t-k >= 2s^2 >= 3/2". Letters real but depending on c: said once at line 32 (regularisation paragraph), echoed in theorem (b), the reading paragraph and gap note.

**T6 / Opus A4, B2; Codex 8.** M3alpha lines 44, 62, 126, 141 and gap (2) (line 234): stale statements replaced by cross-references to Theorem G8B, Corollary G8B', Cor. g8b:cor:inv, Lemma g8b:lem:two; caveats that remain true kept (constant letter [2], positive genus, single angles, depth two suffices, converse criterion, no classical-remainder claim). M2alpha line 31 and gap (i) (line 35): cross-references to Theorem g8c:thm:W3, keeping "not a single iterated integral with constant coefficient". "No referee has read this section" (old M3func line 281, with "and needs a referee") removed; M3func line 3 had no such sentence (the sentence "leaves open item (2)" was replaced by one pointing here). W3int: "it has not been refereed" removed. Nothing about peer review was added.

**T7.** Diffed against main.tex 1500-1580. Changed: Theorem 21 in M4alpha ("21 such terms ... and the constant 1/2"; the "22" remark replaced, with a note that the data file has 22 lines); M3alpha line 8 (R_3: "a rational combination of 74 terms c^e H_w(c) (34 distinct harmonic polylogarithms, of weights two to four)") and the quoted paragraph at line 46 ("Two terms match it ... four components remain. We do not claim that Omega minus these two terms is classical. This is proved at the level of symbols ..."). Prop. 5's displayed statement already matched. Theorem 20 in M2alpha line 15 uses the author-approved sentence, with a note (line 26) that draft 48 has the old wording; same sentence in W3int and the wrapper introduction. main.tex untouched.

**T8 / Codex 1.** Re-run from `proofs/`:
- `g8b_check.py`: pass (exit 0, no False, output identical to stored .out), 1.8 s.
- `g8c_structure.py`: pass (all differences 0), 4.5 s, now with step [7].
- `g8c_numeric.py`: pass (output identical to stored .out), 13.7 s.
- `g8b_pullback.py`: pass, 0.8 s.
Data files now listed with each machine step: M1/M2 none (terms typed in script; g8_zsum.py); M3/M4 g8_classical_dict.json, g8_terms_k2.json (+ g8_cor3_common.py); M5 g8_terms_k3_part1.json, g8_terms_k3_partz2.json; M6/M7 g8_terms_k3_part1.json (run inside scripts/); M8/M9 g8_terms_k4_part1.json, g8_terms_k4_partz2.json; g8b_check.py both k3 files; g8b_pullback.py k3 part1; g8c_structure.py and g8c_numeric.py none. The older scripts M1-M9 were not re-run in this pass.

**T9.** `wrap_g8.tex`: loads the three macro files, order G6 (4 sections), M2alpha, W3int, M3alpha, M3func, M4alpha, introduction updated, `\bibitem{Panzer}` added (CPC 188 (2015) 148-166, arXiv:1403.3385). No duplicate macro conflict (all `\providecommand`; `gcmain`, `gbcor` guarded). `macros_g8b.tex` gained the starred theorem `gbcor`. No `\mathsf`; no math in section headings. Build: 68 pages, 0 Overfull, 0 undefined, 0 "??"; Underfull warnings only. `wrap_g8b.tex` and `wrap_g8c.tex` are superseded (one-line comment added at the top; not rebuilt; they would show undefined cross-references).

## Unsure / changed from the referee
1. Opus A1 step 1 is wrong as written: a small disc around some t_0 in Gamma need not contain infinitely many points of Gamma (Gamma may be a sequence). Repair: disc B around the accumulation point t_* (in U, not a pole of phi), identity theorem for meromorphic functions on B (zeros of a non-zero meromorphic function do not accumulate in B, not even at a pole), then remove the finite set Sigma_t. The statement therefore requires t_* in U, not a pole of phi, phi(t_*) in (0,1), phi(Gamma) in (0,1).
2. Step 5: injectivity of phi_* is proved by G(phi) = const with G(X) = X^{n0}(X-1)^{n1}(X+1)^{n-1} and transcendence of phi, instead of the referee's zero-by-zero argument.
3. The corollary concerns functions C in the class of Definition g8b:def:classical on a simply connected neighbourhood of the arc; real parts and conjugates are not discussed.
4. Genus one of C(e^{i alpha/2}, sqrt(cos(alpha/2))) (example of what is not covered) is checked by hand (z^2 - 2y^2 z + 1 = 0 gives u^2 = y^4 - 1), not by script.
5. Corollary G8B' and the W3int counts are new text and have not been through a second referee.
