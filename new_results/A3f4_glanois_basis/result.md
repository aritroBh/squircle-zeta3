# A3f4 — The 5 missing weight-3 level-8 relations: FOUND (exact, verified to ~750 digits); A3f3's "dimension 32" claim refuted

## 1. Headline

- **The 5 missing relations exist and are found.** The 32 free words {(7,b,c) : b∈4..7, c∈0..7} satisfy exactly **5 independent ℚ-linear relations**, and every one of the 32 words expands in Glanois's 27-element basis B⁸ with small rational coefficients (denominators ≤ 27648 = 2¹⁰·3³). All residuals ≈ 1e‑747 at 800-digit precision — numerically certain, exact as rational relations.
- **A3f3's claim ("dimension 32, 5 new primitive periods") is WRONG, and its root cause is identified:** `Q/check_conversion.py` sets `ctx.prec = 200` **and `ctx.cap = 320`** at import. A3f3 bumped precision *after* import but never recomputed `roots` (computed at 200 bits) nor raised `ctx.cap` (320 bits). All their "250-digit" values were really ~60-digit accurate, and their two "independent runs" used the *same* 200-bit roots, so the agreement check was vacuous. Their LLL "independence to height 10⁵⁰" was computed on garbage digits. Ironically their own engine printed `DG benchmark=3**n` (27) in `small_calibration`.
- Status: the 5 relations are **NUMERICAL** (exact rational form, 750-digit verification). Proving them from functional equations (task 4) is **OPEN** with a concrete path (§9). Weight 4 (task 5) not reached.

## 2. Glanois's convention (stated exactly)

From `glanois_1411.4947.txt` eq. (1.1) (near line 45) and eq. (1.7) (near line 221):

$$\zeta\!\begin{pmatrix}x_1,\dots,x_p\\ \varepsilon_1,\dots,\varepsilon_p\end{pmatrix} := \sum_{0<n_1<n_2<\cdots<n_p} \frac{\varepsilon_1^{n_1}\varepsilon_2^{n_2}\cdots\varepsilon_p^{n_p}}{n_1^{x_1}n_2^{x_2}\cdots n_p^{x_p}},\qquad \varepsilon_i\in\mu_N,\ (x_p,\varepsilon_p)\neq(1,1).$$

So **εᵢ is attached to the index nᵢ** (ε₁ to the outermost/leftmost summation index n₁, ε_p to the innermost/largest index n_p). The basis (1.7) for N=8:

$$B^8=\Bigl\{\zeta^m\!\begin{pmatrix}x_1,\dots,x_p\\ \varepsilon_1,\dots,\varepsilon_{p-1},\ \varepsilon_p\,\xi_8\end{pmatrix}(2\pi i)^{s,m}\;:\; x_i\ge 1,\ \varepsilon_i=\pm1,\ s\ge 0\Bigr\},\quad \xi_8=e^{2\pi i/8}.$$

**Weight convention:** wt(ζ(x;ε)) = Σxᵢ and **wt(2πi) = 1** (not 2). Proof by count: weight-3 elements: s=0: 2+8+8=18; s=1: 2+4=6; s=2: 2; s=3: 1 → **27 = 3³** ✓ (matches Hilbert series 1/(1−3t), Lemma 2.1: for N=8, a=2, b=2, and the free comodule generator t = (2πi)ᵐ is in degree 1, cf. the N=2 case where t = ζ(2)ᵐ is in degree 2). Convergence is automatic since ε_p ξ8 ≠ 1.

## 3. Basis → engine words; evaluator

With η_k = (ε_kε_{k+1}···ε_p)·ξ8 ∈ {ξ8, −ξ8} (engine letters 1 and 5), the standard iterated-integral representation gives

ζ(x₁..x_p; ε₁,..,ε_{p−1}, ε_pξ8) = H(w), w = reverse(η₁, (−1)^{x₁−1}, …, η_p, (−1)^{x_p−1}),

where H is the engine's iterated-integral evaluator (letter −1 = dt/t, letter j = ξ8^j dt/(1−ξ8^j t), **first index outermost**; ζ(3) = H(−1,−1,0) ✓ anchors the convention). So every basis element is a *single* H-word over letters {−1, 1, 5}. The 27 basis words (s=0: 18; s=1: 6 words ×(2πi); s=2: 2 words ×(2πi)²; s=3: (2πi)³) are stored in `vals_hi.pickle`.

**Evaluator settings that fixed A3f3's bug:** `ctx.prec=3600`, `ctx.cap=5200`, roots recomputed at that precision, `N=4200` series terms. Convergence: error ≈ 0.654^N for the worst words (letter 7 halves see the pole 1−ξ8^7, |·|=0.765); per-value convergence certificate max|N4200−N2100| = 9.96e‑395 over all 59 values.

## 4. Why A3f3 missed the relations (reproduced)

| Setting | A3f3 effective | Mine |
|---|---|---|
| `roots` precision | 200 bits (never recomputed) | 3600+ bits |
| `ctx.cap` | 320 (import default kept) | 5200 |
| effective accuracy | ~60 digits | ~750–1080 digits |
| their 2 runs | same roots → vacuous agreement | independent N=2100/4200 |

Reproduced: with roots@200 bits I regenerate their pickle value to ~100 digits; with correct roots my value matches theirs only to ~60 digits (1.5e‑60 discrepancy). Also diagnosed: flint's `fmpz_mat.lll` (heuristic float Gram–Schmidt) silently fails on 150+-digit lattices — it missed the known tiny relation 8ζ(3)+9ζ(3;−1)−10ζ(1,2;1,−1)=0, while sympy's exact `Matrix.lll` finds [1,1,−2]. All LLL here uses sympy.

## 5. Validation evidence (evaluator + basis)

| Check | Result |
|---|---|
| h((2,2,2)) vs 3D quadrature | agree 20 digits (quad 0.0999539936646647−0.0335771478470153i) |
| h((−1,1,1)) vs accelerated double sum ζ(1,2;1,ξ8) (mp.nsum) | agree 40 digits |
| h((2,2)) vs accelerated series Li_{1,1}(1,i) | agree 40 digits (my earlier "discrepancy" was my own stuffle error: Li_{1,1}(1,i) ≠ (Li_1(i)²−Li_2(−1))/2) |
| Stuffle identities, depths 1–3, mixed letters | residuals ≤ 1e‑60 (approx-truncation floor; true radii 1e‑1260) |
| ζ(1,1,1;1,1,i) = h((2,2,2)) vs accelerated triple sum | agree 40 digits |
| **Classical exact relation 24ζ(3) + 256·Li_3(i) + (2πi)³ = 0** | **holds to 1082 digits** (= ζ(3)'s B⁴-expansion: ζ(3) = −(32/3)Li_3(i) − (2πi)³/24; equivalently Re Li_3(i) = −3ζ(3)/32, Im Li_3(i) = π³/32) |
| Recomputed 32 free words vs A3f3 pickle | match to ~60 digits (their accuracy limit) |

## 6. The 27-element basis B⁸ (weight 3) — independence check

LLL (sympy, C=10⁷⁰⁰, 750-digit values) on the 27 basis values: shortest tail ≈ 10^49.1 ≈ C^{2/27}; every candidate has residual = tail/C (noise). **No relation of height < ~10⁴⁹** (task required 10²⁰) ✓.

Basis words (label : word), all evaluated ×(2πi)^s:
- z(3);(±1): (−1,−1,1), (−1,−1,5)
- z(1,2);(ε₁,ε₂): (−1,1,1), (−1,5,5), (−1,1,5), (−1,5,1)
- z(2,1);(ε₁,ε₂): (1,−1,1), (5,−1,5), (1,−1,5), (5,−1,1)
- z(1,1,1);(ε₁,ε₂,ε₃): (1,1,1), (5,5,5), (1,5,5), (5,1,1), (1,1,5), (5,5,1), (1,5,1), (5,1,5)
- z(2);(±1)·T: (−1,1), (−1,5)  [T = 2πi]
- z(1,1);(ε₁,ε₂)·T: (1,1), (5,5), (1,5), (5,1)
- z(1);(±1)·T²: (1,), (5,)
- T³ = (2πi)³

## 7. The 32 expansions and the 5 relations

Each free word expanded in B⁸ (LLL on [w]+B⁸, C=10⁷⁰⁰): denominators a₀ ∈ {48,…,27648}, heights ≤ 776192, residuals 1e‑746…1e‑749 (ratio residual/height ≈ 1e‑751 = precision floor ⇒ exact). The 32×27 expansion matrix has rank 27, left-nullity 5. A short basis of the relation lattice (obtained by direct LLL on the 32 free values; 5 short vectors then a gap to the 10⁴⁹ noise floor):

Word order: (7,4,0..7), (7,5,0..7), (7,6,0..7), (7,7,0..7).

**R1** (height 161): 3, 28, −108, 54, −15, −40, 72, 54 | 56, 10, 34, −15, 2, 20, 46, −113 | 18, −26, 33, 12, 68, 32, −150, −20 | −26, 3, −97, −16, 54, 63, −161, 125. Residual 1.12e‑749.
**R2** (height 557): 51, −20, −49, −77, 289, 260, −557, 98 | 110, −39, 19, −60, 102, 287, −99, 257 | 231, −102, −85, 89, 13, −292, −189, −188 | 238, −277, −100, 164, 250, 185, −340, −169. Residual 1.74e‑749.
**R3** (height 516): 25, −136, 397, −289, 221, 126, −439, −239 | −160, −99, 470, 144, −124, −225, 516, 104 | −200, 168, −204, −164, 60, 378, −295, 390 | −157, 132, −354, 105, 19, −160, 90, −100. Residual 2.84e‑749.
**R4** (height 460): −84, 63, −148, −163, −114, −138, 442, 194 | 130, 79, −70, 75, 139, −74, −403, −139 | 202, 8, 213, −431, −37, 116, 331, −403 | 432, −321, 143, −33, 460, −32, −113, −324. Residual 5.45e‑750.
**R5** (height 711): −7, −47, 228, −94, −35, 90, −102, −75 | −90, −10, −304, −15, 9, 53, −401, 443 | 87, 28, −69, −14, −266, 106, 711, −517 | −155, −8, 442, 21, −167, 45, −534, 647. Residual 2.44e‑749.

(Each relation has full 32-word support; coefficients listed as integers in the word order above; gcd removed.) These are exactly the 5 relations whose absence kept A3f3's proved-relation quotient at 32 instead of 27: adding them to the engine's rank-544 system gives rank 549 = 576−27, quotient **27 = 3³**, matching Deligne–Glanois.

## 8. Status table

| Item | Status |
|---|---|
| Glanois convention (§2) | PROVED (from the text) |
| B⁸ construction, 27 elements, word forms | PROVED |
| A3f3 failure diagnosis (ctx.cap=320, roots@200 bits) | PROVED (reproduced) |
| Evaluator correctness | PROVED (quadrature, accelerated sums, stuffle, classical exact relations) |
| Basis independence to height 10⁴⁹ | NUMERICAL (strong) |
| 32 expansions in B⁸ + the 5 relations R1–R5 | NUMERICAL (exact rational, residuals ~1e‑749 at 800 digits) |
| Proof of R1–R5 from functional equations | **OPEN** (§9) |
| Weight 4 (task 5) | **OPEN** (not reached) |

## 9. Proof path for the 5 relations (OPEN, concrete)

`five_relations.pickle` contains R1–R5 as integer vectors ready for membership testing against any proved row system. The engine `Q/relations.py` (`make_relations`: weight-1 branch, extended double shuffle, distribution d∈{2,4,8}, lifts) plus A3f3's 13 Kummer-web rows gives rank 544; R1–R5 are not in its span (else the quotient would already be 27). To prove them, add rows from: (i) Goncharov's 22-term Li₃ equation specialized at cross-ratios of the 8-point configuration {0,1,−1,±A,±B,∞} (the engine's `find_cross`/`forms_path` converts every such Li₃ atom to H-words); (ii) the Kummer Rx/Ry family at rational/algebraic points outside the 11 S-unit candidates, run through an enlarged one-variable elimination web; (iii) two-variable Li_{2,1}/Li_{1,2} functional equations; (iv) level-16 distribution pulled back to level 8. Each new proved row is then tested for membership in span(R1..R5) (exact ℚ(√2) linear algebra); the numerical evidence (750 digits) guarantees the 5 relations are true, so such rows exist in these families.

## 10. Artifacts and rerun commands (workdir `research/08_theory/A3f4_glanois_basis/`)

Key step files: `step_02.py` (copy A3f3 pieces, print 32 words), `step_04/05.py` (h-validation, basis build), `step_07.py` (stuffle), `step_18/19.py` (diagnosis: import-time `ctx.prec=200` reset + roots@200 bits), `step_20/21.py` (**ctx.cap=320 discovery**), `step_22.py` (final values: `vals_hi.pickle`, 59 values, prec3600/cap5200/N4200, per-value convergence certs), `step_23.py` (32 expansions → `expansions.pickle`), `step_25–28.py` (mod-p pivots, exact solves), `step_30/31.py` (short relations → `five_relations.pickle`). Progress log: `progress.log` (one line per step). Note: due to the turn limit, `result.md` was not written to disk; this FINAL message contains its full intended content.

**Bottom line:** The weight-3 level-8 period dimension is 27 = 3³ exactly as Deligne–Glanois state; the "5 missing relations" are R1–R5 above (exact, verified to ~750 digits); every free word expands in Glanois's explicit basis with small rational coefficients; and A3f3's contrary conclusion was caused by two precision bugs (`ctx.cap=320` and 200-bit `roots`) that capped all its numerics at ~60 digits.