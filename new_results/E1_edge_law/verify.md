# E1 verification: adversarial check of the edge law

Verifier's scripts and outputs: `verify/` in this folder. They were run from a copy in
`<scratch>/verify_E1/`; the saved copies differ only in their paths. Every number below is printed
in a log under `verify/out/` or `verify/rerun_logs/`.
Labels: **PROVED** means exact symbolic computation or a checked written argument. **NUMERICAL** means agreement at high
precision, residual printed. **CONJECTURE** means a pattern with no proof.

## Verdict: PARTIALLY CONFIRMED

Confirmed:
- the formula for Z_e, read as the ζ(3) coefficient of the *local* part of each edge integrand;
- the radial-projection method;
- every reduction (C_m, rhombus, parallelogram, affine-regular hexagon, tangential and angle forms);
- the exact closed forms for fully local polygons;
- the correction to the C_m asymptotics.

**Not confirmed, and false as an absolute statement:** "the ζ(3)-part of P‴(1) is Σ_e Z_e for every centrally
symmetric convex polygon", and "the non-cross terms cancel (PROVED for all m)".

I built a rational hexagon that is *not* fully local. All its constants reduce to alternating MZVs of weight ≤ 3, so
its ζ(3) coefficient is canonical. Even so, P‴(1) contains (159√5/125)ζ(3) more than the edge law predicts. The extra
ζ(3) sits in the non-local "corner log × other-kink log" terms. The solver's §3 theorem says only that no *individual
term* of their series contains ζ(3), and that is true. It does not say that the sum is ζ(3)-free, and here the sum is
not.

## 1. Claim-by-claim status

| solver claim | status after verification |
|---|---|
| Lemma 1 (exact arclength in the radial parametrisation), Lemma 2 (polynomial × logs of affine functions) | **CONFIRMED.** I re-derived it in sympy straight from x_p = y·F^{-1/p} (`t5`). Edge by edge, its ∫E₃ equals my polar-angle jet integral to 40 digits. |
| g₂₁, g₁₂ and Z_e = 12∫(g₂₁+g₁₂) = 6αβ(h²/L)(3(αu+βv)−2) as the ζ(3) coefficient of the local (J = 0) part | **CONFIRMED (PROVED).** My own derivation gives local ζ(3) coefficient 18√5/125 on the r = 2 edge of my hexagon, and the formula gives the same. All the solver's logs re-run identically. |
| Equivalent forms (α = L sin ε₀/(2h₋), angle form, tangential form) and reductions to C_m, 24a²b²/r³, 6h²/L, −(3/4)h²/L | **CONFIRMED.** Re-runs of s10, s14, s16 are identical. I hand-checked the angle form and the rhombus case of the tangential form. |
| Closed forms: parallelograms, affine-regular hexagons, H_p, and D15's hexagon now a theorem | **CONFIRMED.** Exact (s04, s08) and independently NUMERICAL (my jet code agrees to ≤ 4.4e-55; my own affine-regular hexagon agrees to 9e-55; my own fixed-weight parallelogram is identified by PSLQ, §3.2). |
| §3 "non-cross cancellation, PROVED for all m" | **The literal statement is true; the intended meaning is false.** No term of the canonical series contains ζ(3), but the series can sum to a nonzero multiple of ζ(3). §3.3 gives an explicit non-degenerate polygon where this happens. |
| §0.1 "P‴(1)\|_{ζ(3)} = Σ Z_e, every centrally symmetric convex polygon" | **TRUE ONLY** for fully local polygons (parallelograms, affine-regular hexagons), or as a definition of the "canonical" coefficient. **FALSE** as an absolute coefficient in general (§3.3). |
| §5 caveat | **Points the right way, but understated.** The hidden ζ(3) is not a degeneration effect (r → 0) and not only a basis ambiguity (regular m ≥ 4). It is present at r = 2, deep inside the disc of convergence, where the basis is canonical. |
| §8 paper paragraph: "Adding one Z per side gives all the ζ(3) in P‴(1)"; "ζ(3) shows up only in the mixed terms" | **FALSE in general.** Must be restricted to fully local polygons. |
| §7.1 C_m = −3π⁴/m³ + 25π⁶/(8m⁵) + O(m⁻⁷); ledger's 9π⁵/(16m⁴) wrong | **CONFIRMED.** My series (`t8`): C_m = −3π⁴/m³ + 25π⁶/(8m⁵) − 609π⁸/(640m⁷) + …. Only odd powers occur, so O(m⁻⁷) is right. |
| §2.2 Z_e independent of λ | **CONFIRMED** for the canonical coefficient. It also holds for the absolute coefficient whenever the λ_j are logs of algebraic numbers: a constant shift of the ℓ_j only lowers the log-degree. Trivial exception: λ_j that themselves contain ζ(3). |

## 2. Re-runs (copy in `<scratch>/verify_E1/E1_copy`)
`run_all.sh` was re-run from the copy, together with s13 (40 digits), s16, s17 and s06 at 60 digits. All 19 logs are
byte-identical to the solver's: s01 to s17, s06 at 30 and 60 digits, s11b. Copies are in `verify/rerun_logs/`.

## 3. Independent computations (the verifier's own code)
`verify/vjet.py` parametrises the boundary by the polar angle, ρ(θ) = F_p(u_θ)^{-1/p}. It carries the d-dependence by
*generic truncated power-series arithmetic*. It uses neither the solver's E_k nor D15b's hand-expanded series.
`verify/t4_fd.py` uses no d-expansion at all. It takes the direct arc length at p = 1 + kh, then a symmetric
(2K+1)-point polynomial fit.

### 3.1 Sanity (`t1`, 60 digits)
|polygon|quantity|difference from the solver's exact closed form|
|---|---|---|
|regular hexagon|P‴(1)|3.1e-55|
|H_p, power convention|P‴(1)|1.8e-55|
|rhombus a = 1, b = 2|P‴(1)|4.4e-55|

### 3.2 Fully local, non-regular, fixed weights (`t2`): K_p = {(|x|^p + |x+2y|^p)/2 ≤ 1}
Vertices (2,−1), (0,1), (−2,1), (0,−1). The sides are 2√2 and 2, so this is not a rhombus. Per-arc PSLQ of P‴_e/L_e
on {ζ(3), π²log2, log³2, π², log²2, log2, 1} gives the same relation at 60 and at 110 digits:
- arc 1: [16, −6, 1, 16, 2, 96, 96, 0], so the ζ(3) coefficient is 3√2/4 = 6h²/L (difference 1.6e-61 at 60 digits);
- arc 2: [4, −6, 1, 4, 2, 24, 24, 0], so the ζ(3) coefficient is 3 = 6h²/L.

The edge law holds absolutely here, as it should for a fully local polygon.

### 3.3 Not fully local: the rational hexagon (`t3`, `t4`, `t5`, `t6`, `t7`)
**Definition.** K_p = {|x|^p + |y|^p + 3|x−y|^p ≤ 4^p}. At p = 1 the vertices are (1,0), (2,2), (0,1) and their
negatives. Every |⟨v_j, V_k⟩| is 1/4 or 1/2, so only log 2 enters the constants μ.

**The non-local ratios r** are 2, 1/2 and 1 on the three edges. So each edge integral is an alternating MZV of weight
≤ 3 times L_e. In that algebra the only indecomposable of weight 3 is ζ(3). The basis is
{ζ(3), π²log2, log³2, π², log²2, log2, 1}. So the ζ(3) coefficient is canonical, modulo the standard (unproven)
ℚ-linear independence of these constants. At the motivic level it is canonical unconditionally.

| arc (edge) | r | edge law Z_e | ζ(3) coefficient found by PSLQ (60, 110, 150 digits, polar jet; and on the FD-only values at 38 and 60 digits) | relation on [P‴_e/L_e, ζ(3), π²log2, log³2, π², log²2, log2, 1] |
|---|---|---|---|---|
| 1: (1,0)→(2,2) | 2 | 18√5/125 | **231√5/500** | [500, −231, −111, −52, 80, 1440, −2400, 0], residual 2.5e-141 at 150 digits |
| 2: (2,2)→(0,1) | 1/2 | 18√5/125 | **231√5/500** | same (edge 2 is the mirror image of edge 1) |
| 3: (0,1)→(−1,0) | 1 | −15√2/128 | −15√2/128 | [128, 15, −4, −1024, 4, 3072, −1536, 0] |

On arcs 1 and 2, the remainder (P‴_e − Z_e ζ(3))/L_e has **no** relation on the ζ(3)-free basis at either precision.
PSLQ returns only noise with coefficients near 1e7 and then 1e14. On arc 3 the remainder is elementary:
32·R = π²log2 + 256log³2 − π² − 768log²2 + 384log2.

**Total, NUMERICAL.** PSLQ at 110 and at 150 digits gives the same relation, residual 8e-141:

P‴(1) = √5[(231/125)ζ(3) + (111/125)π²log2 + (52/125)log³2 − (16/25)π² − (288/25)log²2 + (96/5)log2]
      + √2[−(15/64)ζ(3) + (1/16)π²log2 + 16log³2 − (1/16)π² − 48log²2 + 24log2]
      = 19.900266363148066980424545474916156477954312358839…

The finite-difference value (direct arc length only, 170 digits, h = 4e-4, K = 26) equals this closed form at every
digit printed. On arc 1 alone it agrees to 1.9e-63. If the closed form used the edge-law coefficient instead, it would
miss by (159√5/125)ζ(3) = 3.419.

The true ζ(3) coefficient is (231/125)√5 − (15/64)√2. The edge law gives (72/125)√5 − (15/64)√2.

**Where the extra ζ(3) lives (`t5`).** I derived E₃ symbolically on edge 1 from x_p = y·F^{-1/p}. Its local part
(PROVED, exact Beta moments via sympy derivatives of Γ) has ζ(3) coefficient exactly 18√5/125 = Z_e. The edge has
40 non-local monomials t^n log^a t log^b(1−t) log^c(1+t), c ≥ 1. Each was identified by PSLQ (NUMERICAL), and
together they carry **159√5/500 ζ(3)**. The sum is 231√5/500, which matches the per-arc PSLQ. Their sum, together with
the local part, also reproduces the polar arc value to 40 digits.

**Robustness.**
- *Task's fixed-weight convention* (`t7`): K_p = {(1/4)|y|^p + (3√2/4)|(y−x)/√2|^p + (1/4)|x|^p ≤ 1}, the same K_1.
  It gives the same per-arc ζ(3) coefficients (231√5/500, and −15√2/128 on arc 3) at 70 and at 110 digits, although
  P‴ itself differs (total −85.309…).
- *Sheared hexagon* (b = 1; vertices (1,0), (4,2), (1,1)), which has no mirror symmetry:
  - arc 1: PSLQ gives 615√13/8788, the edge law 99√13/4394;
  - arc 2: PSLQ gives 267√10/2000, the edge law 99√10/2000;
  - arc 3 (r = 1) matches: −3√5/160.

  The same relations hold at 60 and 110 digits (residuals 1e-51 and 1e-100).

So the "edge law" is exactly the ζ(3) coefficient of the local terms. In general it is not the ζ(3) coefficient of
P‴(1).

## 4. Reading of the PROVED steps
- **Lemma 1, the smoothness argument, Lemma 2 and the Beta lemma:** sound. The bound t^{−δ}(1−t)^{−δ}·(logs)^k
  justifies differentiating under the integral sign. The Beta table matches my independent Γ-derivative evaluation.
- **§3 Theorem:** correct as stated. Each term of A_e, B_e lies in R_e[log r_i]·(Beta moments with a+b ≤ 2), and
  those moments lie in ℚ + ℚζ(2). But this says nothing about the ζ(3) content of the *sums* A_e, B_e. The solver's
  own s12 example, ∫log²t·log(1+t) = 4log2 − 6 + π²/6 + (3/2)ζ(3), is precisely a monomial that occurs on an r = 2
  edge. The gap is in the inference from the theorem to "ζ(3) comes only from the cross terms" (§0.2, §3 heading,
  §8): it does not follow, and §3.3 shows it is false.
- **§3's answer to the resummation caveat** (local coefficients are polynomials, so the Beta rule is exact) is right
  for the local part only. It does not cover the non-local part, which is where the hidden ζ(3) sits.
- **§5:** "the hidden ζ(3) emerges as r → 0, where r − 1 reaches the edge of its disc" is misleading. At r = 2 the
  expansion variable is 1/r − 1 = −1/2, and the hidden ζ(3) is already there, at a non-degenerate polygon.
- **§7.2** ("PROVED … in the canonical sense for all m, and as an absolute statement for m = 2, 3 regular and all
  fully local polygons") is the accurate formulation. The headline (§0.1), the summary claim 3 and §8 should be
  brought into line with it.
- **Regular 2m-gons, m ≥ 4:** C_m is the canonical (local) coefficient. In the light of §3.3, whether P‴(1) has
  ζ(3) coefficient C_m in any natural basis is **open**, and the answer is probably no. This is more than a choice of
  basis, because the non-local terms generically carry their own ζ(3)-type indecomposables. D15b/D15c's "extraction
  theorem" for m ≥ 4 is in the same position.

## 5. Corrections to make
1. In result.md §0.1, the summary and §8, replace "the ζ(3)-part of P‴(1) is Σ_e Z_e (every polygon)" with: "Σ_e Z_e is
   the ζ(3) coefficient of the local (corner–corner) terms. It equals the ζ(3) coefficient of P‴(1) for
   parallelograms and affine-regular hexagons (PROVED). In general it does not: for {|x|^p+|y|^p+3|x−y|^p ≤ 4^p} the
   coefficient is (231/125)√5 − (15/64)√2, against the edge law's (72/125)√5 − (15/64)√2 (NUMERICAL, PSLQ at 110 and
   150 digits, plus finite differences)."
2. Rename "non-cross cancellation" to "termwise absence of ζ(3) in the canonical expansion". Keep the true and useful
   part: single-corner log³ terms are rational, log² terms give ℚ + ℚζ(2), and the polar-angle ζ(3) is an artefact.
3. In the ledger, D15c's "cross = total" stays PROVED only for fully local polygons. For m ≥ 4 regular polygons it is
   a statement about the canonical coefficient.
4. The C_m asymptotics correction stands: C_m = −3π⁴/m³ + 25π⁶/(8m⁵) + O(m⁻⁷).

## 6. Files (`verify/`)
- `vjet.py`: jet perimeter in the polar angle and the direct perimeter.
- `t1_sanity.py`: against the closed forms.
- `t2_parallelogram.py`: fixed-weight parallelogram PSLQ.
- `t3_hexagon.py`: rational hexagon, per-arc PSLQ. The argument b is the shear: 0 or 1.
- `t4_fd.py`: finite differences of the direct arc length.
- `t5_decomp.py`: independent symbolic E₃, local versus non-local ζ(3).
- `t6_pslq_fd.py`: PSLQ on the finite-difference values.
- `t7_fixed_and_total.py`, `t7b_total150.py`: task convention and total-perimeter PSLQ.
- `t8_misc.py`: finite-difference total against the closed form, the affine-regular hexagon law, the C_m series.
- `out/`: all logs and JSON.
- `rerun_logs/`: the solver's logs, re-run.
