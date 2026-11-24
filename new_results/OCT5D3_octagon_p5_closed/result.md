OCT5D3 — the closed form of the regular octagon's P⁽⁵⁾(1)

# OCT5D3 — the closed form of the regular octagon's order‑five coefficient

**Bottom line.** All four tasks are complete. P⁽⁵⁾(1) of the regular octagon is now in closed form: P⁽⁵⁾(1)/u = Wt1+Wt2+Wt3+Wt4+Wt5 (Wt0 = 0 exactly, u = √(2+√2)), with Wt1–Wt4 in classical constants, and Wt5 = (52 classical M4f atoms) + (9 Li₃,₂) + (21 invisible products). The weight‑5 symbol identity and obstruction are **PROVED** (exact symbol algebra); the value‑level identity is **NUMERICAL**, verified at two precisions (1100/1500 digits) with residuals at working precision and a PROVED symbol check of the whole closed form. It reproduces p5_800's P⁽⁵⁾(1) to 4.9e‑794. `verify_oct5d3.py` runs in **52.6 s** at 360 digits: **25 PASS lines, rc 0**.

Two findings of independent value:
1. **OCT5D2's stored Li₃,₂ values (`atom_vals430.pkl`) are off by a factor 1/x** (off‑by‑one series) — a bug that sank OCT5D2's closure attempt. The value consistent with the symbol convention gsym(0,0,1/x,0,1/(xy);1) is the genuine double polylog Li₃,₂(x,y) = Σ_{n>m>0} xⁿyᵐ/(n³m²): Li₃,₂(b,1) = 0.02554867181188595195138362416354460326083… (= b×their stored value to 3.34e‑423).
2. **The "weight‑5 part" is motivic (top G‑words), not monomial‑degree.** The degree‑5 block (−7.375079049848819…, matches OCT5D2's w5_430 to 4.84e‑420) is mixed‑weight; the pure weight‑5 part is **Wt5 = 76.230816325968618891401147030033662…**, whose symbol is the fitted one.

## Task 1 — classical representative in the M4f basis (PROVED, exact)

M4f basis: 70 atoms (14 Li₅ at {b, −b, −b², 1/√2, 1/2, 1−1/√2, 2−√2} ∪ extras {2b, (2−√2)/4, (2+√2)/4, −b/2, (4−3√2)/8, b⁴, 1−b⁴}; 54 products; ζ(5), L(χ₈,5)); 52 have nonzero symbol. (One correction to my first pass: M4f's Li₄ argument is −b², matching R5_best.json and the paper's eq:R.)

- C = (−3600+2160√2)Li₃,₂(b,1) + (11880−8280√2)Li₃,₂(1/√2,1): its (3,2) pair = **(−1485+1035√2, 900−540√2)** = the obstruction (exact).
- T = S(W₅ᵒᶜᵗ) − S(C): obstruction (0,0); lies **exactly** in the span of the 52 visible atoms: rank 30/32 over ℚ(√2), 0 inconsistent rows, **residual tensor exactly 0**, kernel dim 22. RREF particular solution (30 terms, e.g. Li5(b): −17523420/23+12318650/23√2; full list in result.md / verify output).

## Task 2 — precision (NUMERICAL)

New parallel engine `gengine2p.py` (same mathematics as OCT5D2's validated gengine: exact‑form series at t = 1/2 + Chen concatenation; one worker per length‑3 subtree). **Time per level (printed):** (dps 90, N 320) 7.8 s serial; (430, 1440) 66.1 s serial / 9.0 s parallel, residual 4.21e‑432 (matches old table); (1450, 560) 83.7 s serial; **(1650, 5700) 184.5 s parallel, residual 3.54e‑1652**. Atoms (14 Li₅, 7 Li₄, 3 Li₃, Li₂(b), ζ(2..5), L(χ₈,3), L(χ₈,5), 9 Li₃,₂) at 1650 digits, cross‑validated (duplication 2.8e‑1629; L(χ₈,5) vs M4e's certified value 7.09e‑1652). True weight parts (term weight = |G‑word| + e + f), dps 1650:
```
Wt0 = 0 (EXACT, QS2)   Wt1 = 178.61011311121766228088213936448409...
Wt2 = -532.49879617214258393968401364283616...  Wt3 = 314.42432505669915486229105783717946...
Wt4 = 52.058844836112780766819569870462025...  Wt5 = 76.230816325968618891401147030033662...
sum = 88.825303157855632861709900459323071...  (vs p5_800 P5/u: 8.24e-800)
```

## Lower weights (Wt1 PROVED; Wt2/3/4 NUMERICAL, identical at two precisions)

- **Wt1 = (4800−3360√2)λ + (2880−1920√2)L** (derived exactly; LLL residual 0.0)
- **Wt2 = (360−280√2)π² + (−20880+14760√2)λ² + (−19200+13440√2)λL + (−7680+5280√2)L² + (−3840+2880√2)Li₂(b)** (identical 800/1500, h=2.09e4)
- **Wt3** (12‑atom paper basis): (−4320+2880√2)Li₃(b) + (15120−10800√2)Li₃(−b) + (11520−8640√2)Li₃(1/√2) + (−17985/2+6915√2)ζ(3) + (6480−4560√2)L³ + (17280−12240√2)λL² + (34200−24120√2)λ²L + (34410−24330√2)λ³ + (5760−4320√2)Li₂(b)λ + (5760−4320√2)Li₂(b)L + (1370−980√2)π²λ + (−900+660√2)π²L (identical 900/1500, h=6.88e4)
- **Wt4** (29 atoms, paper R‑basis; coefficients in coeff_tables.py / closed_form_section.md; identical 900/1500, h=1.37e6)

## Task 3 — the invisible closure (NUMERICAL, accepted)

Basis (51 atoms): 21 planner invisible (ζ(2)×12 weight‑3, ζ(3)×5, ζ(4)×2, ζ(5), L(χ₈,5)) + **K32** (= Li₃,₂(−b²,1) − 8Li₃,₂(b,1) − 32Li₃,₂(−b,1) − 16Li₃,₂(b,−1) − 16Li₃,₂(−b,−1); the Li₃,₂ symbol kernel is PROVED dim 1) + **7 mixed Li₃,₂×classical kernel constants** (PROVED: the (52+9)‑atom symbol kernel is 29 = 22 classical + K32 + 6–7 mixed) + 22 classical symbol‑kernel combos.

- LLL at **1100 digits**: rel1100, height 2.13e9, residual 1.26e‑1090 (working precision).
- LLL at **1500 digits**: rel1500, height 2.2e9, residual 1.47e‑1490 (working precision).
- rel1500 **verified by direct evaluation** at 900 (1.17e‑894), 1100, 1500 (1.4e‑1494).
- rel1500 − rel1100 lies in the **verified 26‑dim atom value‑kernel** (exact rank check): the two runs coincide modulo exact atom‑relations. **ACCEPTED.**
- Height: lcm denominator **8832 = 2⁷·3·23**, max |numerator| **281,424,457** (LLL integer height 2.2e9 with the target coefficient). **L(χ₈,5) coefficient: −400/3 + (39752/3)√2** (basis‑dependent).
- **Final PROVED check: S(closed form of Wt5) = S(W₅ᵒᶜᵗ) exactly**, obstruction = (−1485+1035√2, 900−540√2).

## Task 4 — the whole number + local law

P⁽⁵⁾(1) = u·ΣWt_k = **164.1277591133059706510456733124543613051…** — matches p5_800 to 4.88e‑794 (P5/u to 2.64e‑794). prop:order5 (m=4, formulas cross‑checked vs E5_edge_law):
```
C_4^(5)   = (33750-23805√2) cos(π/8)  = 78.2028433737847736987918987074
C_4^(2,3) = (-17370+12315√2) cos(π/8) = 42.5354327315259866080186809336
closed form (basis-dependent):  ζ(5):        u·(82937921/2944 - 56662239/2944√2)  = 1760.7776471224051...
                                ζ(2)ζ(3):    u·(28661199/2944 + 4602059/2944√2)   = 22073.6243104001168... (π²ζ₃ folded in)
non-local parts (basis-dependent): ζ(5): 1682.5748037486...,  ζ(2)ζ(3): 22031.0888776686...
```
The octagon is far from fully local (consistent with the nonzero depth‑two obstruction); these non‑local parts shift under the 26‑dim atom value‑kernel — **basis‑dependent, as stated**.

## The closed form (full table in `closed_form_section.md` / `coeff_tables.py`; λ = log 2, L = log(1+√2), b = √2−1)

Wt1–Wt4 as above. **Wt5** = 52 classical atoms (coefficients β incl. all 14 Li₅'s, e.g. Li5(b): 80771059/138−21711787/92√2; Li5(b⁴): 13430191/1472−64746761/8832√2; λ⁵: 4902637/192−3467879/192√2; …) + 9 Li₃,₂ (e.g. Li32(b,1): 1737787/23−2530311/46√2; Li32(−b²,1): −16995/4+13255/4√2) + 21 invisible products (ζ(5): 82937921/2944−56662239/2944√2; L(χ₈,5): −400/3+39752/3√2; ζ(2)ζ(3): −16159155/2944−8342511/2944√2; …). Complete coefficients: `coeff_tables.py` (145 lines), machine‑generated.

## Evidence table (all numbers printed from code)

| # | Check | Result | Residual | Digits | Producer |
|---|---|---|---|---:|---|
| 1 | obstruction S(W₅ᵒᶜᵗ) | PROVED exact | 0 | — | verify_oct5d3.py §A |
| 2 | C's pair = obstruction | PROVED exact | 0 | — | verify §A |
| 3 | S(W₅)−S(C) in 52 M4f atoms span, residual 0 | PROVED exact | rank 30, kernel 22 | — | verify §A |
| 4 | Li₃,₂ symbol kernel dim 1 (K32) | PROVED exact | 0 | — | verify §A |
| 5 | S(closed form Wt5) = S(W₅ᵒᶜᵗ) | PROVED exact | 0, obstruction matches | — | wt5_primitive step |
| 6 | G‑table dps 1650 N 5700 | built 184.5 s | 3.54e‑1652 | 1652 | gengine2p, build_big.log |
| 7 | Wt0 = 0 | PROVED (QS2) | 0 | — | verify §B |
| 8 | ΣWt_k vs p5_800 P5/u | NUMERICAL | 8.24e‑800 | 800 | assembly step |
| 9 | Wt1 exact | PROVED + num | 0.0 (dps 800) | 800 | Wt1 step; verify §B |
| 10 | Wt2, Wt3, Wt4 identical 800(900)/1500 | NUMERICAL | 1.4e‑1497 / 1.9e‑1497 / 1.5e‑1495 | 1497 | q2_lll steps; verify §B |
| 11 | Wt5 closure, direct eval 900/1100/1500 | NUMERICAL | 1.17e‑894 / 1.26e‑1090 / 1.4e‑1494 | 1490 | closure steps; verify §B (1.1e‑351) |
| 12 | rel1500 − rel1100 ∈ atom‑kernel | exact rank 26=26 | — | — | membership step |
| 13 | closed form reproduces P⁽⁵⁾(1) | NUMERICAL | 4.88e‑794 | 794 | assembly; verify §C |
| 14 | Li₃,₂ convention (OCT5D2's 1/x bug) | NUMERICAL | 3.34e‑423 | 423 | atoms step; verify §D |
| 15 | duplication Li₅(b⁴) | NUMERICAL | 2.81e‑1629 | 1629 | atoms step; verify §D |
| 16 | L(χ₈,5) vs M4e | NUMERICAL | 7.09e‑1652 | 1652 | atoms step; verify §D |
| 17 | verify_oct5d3.py end‑to‑end | 25 PASS, rc 0, 52.6 s | 1e‑349…1e‑363 | 360 | verify run |

## Files in `new_results/OCT5D3_octagon_p5_closed/`

- **`verify_oct5d3.py`** — deliverable verifier (52.6 s fresh at 360 digits, 25 PASS, rc 0; reproduces symbol algebra from parsed_mons, builds the G‑table, assembles all weight parts, checks every closed form, prints the local‑law comparison).
- Modules: `oct5d3_basis.py` (M4f basis + C + value specs), `oct5d3_atoms.py` (atom evaluator, corrected Li₃,₂, K32), `gengine2.py`/`gengine2p.py` (serial/parallel G‑table), `q2_lll.py` (ℚ(√2) LLL), `coeff_tables.py` (all closed‑form coefficients), `closed_form_section.md`, `result.md` (full report with the complete closed form).
- Data: `gtab2_d1650_N5700.pkl` (G‑table ~1600 digits), `wparts_true_1650.pkl`, `atoms_1650.pkl`, `symfit_m4f.pkl`, `kernel_red.pkl`, `kernel61.pkl`, `mixed61_vals.pkl`, `atom_kernel.pkl` (26 verified atom relations), `rel1500_frac.pkl`/`rel1100_frac.pkl`, `wt5_primitive.pkl`, `closed_forms_lower.pkl`, `rewrite_all.pkl`, `build_big.log`, and the auto‑saved `step_01.py…step_NN.py` + `transcript.json` (full computation transcript; key producers: symbol fit step (symfit_m4f), gengine2p build step (build_big), assembly/Wt steps, closure steps (rel1100/rel1500), kernel‑membership step, verify run).

**Labels recap.** PROVED (exact symbol/QS2 algebra): obstruction, C's pair, the exact identity S(W₅) = S(classical_M4f) + S(C) with the printed representative, the K32 and 29‑dim kernel structures, Wt0 = 0, Wt1, and S(closed form Wt5) = S(W₅ᵒᶜᵗ). NUMERICAL (digits, residuals printed): the 1650‑digit G‑table and weight parts, Wt2–Wt4 two‑precision forms, the Wt5 invisible closure (two‑precision + kernel‑membership acceptance), the P⁽⁵⁾(1) assembly check, and the atom cross‑checks. The exactness of the value‑level identity Wt5 = (closed form) is NUMERICAL (modulo the period conjecture), as is standard for this kind of result; the symbol‑level content is exact.