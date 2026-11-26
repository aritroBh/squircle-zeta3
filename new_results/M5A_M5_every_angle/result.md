# M5A: M₅ at every angle

**Summary.** M₅(α) is now known at every angle in the same form as M₄. Every zeta part is a rational combination of terms c^e H_w(±c) with −4 ≤ e ≤ 1, and each part is **PROVED for every Taylor coefficient** by the same Z-sum method as G2/G8. The new parts that the paper asks for are [ζ(2)]M₅ (60 terms and the constant ½) and [ζ(2)²]M₅, which has a short closed form. All four ladder predictions agree with independent fits. The checks against the definition hold to about 52 digits at six angles, and to 10⁻¹⁷⁰ at the square corner.

Work dir: `new_results/M5A_M5_every_angle/`. Python: project `.venv`. Labels: **PROVED** (exact algebra, all n), **EXACT THROUGH c^N**, **NUMERICAL** (digits). Every number below was printed by a script in this folder (output files `*.out`). No existing file was modified. The library and data files that I used were copied here: `g8_allN.py`, `g8_zsum.py`, `g6_common.py`, `g8_terms_k3_partz2.json` (Ψ), `g8_terms_k4_partz2.json` (Ψ₂) and `g8_terms_k4_part1.json` (R₄, only for the fit self-test).

## 1. Notation

c = cos(α/2), and H_w are harmonic polylogarithms with letters 0, ±1, as in the paper. "c^e H_w(±c)" means the pair c^e(H_w(c) + s H_{−w}(c)) with s = (−1)^{e+ν(w)}. This equals c^e H_w(c) + (−c)^e H_w(−c) and is even in c. The JSON lines are `["c^e He[w]", q]`, as in G8, plus the constant lines `"1"` and `"c^-2"`. Shorthands:

- A = artanh c, L = log(1 − c²);
- D(c) = H_{0,1}(c) + H_{0,−1}(c) − H_{−1,1}(c) + H_{1,−1}(c) = Li₂(c) − Li₂(−c) − H_{−1,1}(c) + H_{1,−1}(c);
- Φ = ½ − (1−c²)/(2c)·A (paper, Thm flat), F₂ = [ζ(2)²]M₄ (G8), Ψ = [ζ(2)]M₃, Ψ₂ = [ζ(2)]M₄.

Throughout, [z] means the formal coefficient in the sense of Definition g6:def:formal: each ψ^(m−1) at an integer contributes ζ(m) as a symbol, and products come from products. This is the convention of Theorem thm:zladder, whose kernel is log Γ(1+ξ) + log Γ(1+η) − log Γ(1+u) = Σ (−1)^m ζ(m)/m·(ξ^m + η^m − u^m). So ζ(4) and ζ(2)² are separate monomials.

The convention matters for the ladder rung [ζ4]M₅ = Ψ + G, for the values of G and of Ψ + G taken separately, and for the statement that [ζ(4)]M₅ contains D. It does not matter for the π⁴ bracket π⁴(Ψ/90 + 7G/180), and it does not matter for the ζ(5) and ζ(2)ζ(3) parts. The referee's own derivation from Beta functions at integer points gives exactly the same split for every n ≤ 150. So the split is the same on these two natural routes, but it is still a convention.

## 2. The theorem

**Theorem (M₅ at every angle; PROVED).** For 0 < α < π,

M₅(α) = R₅(c) + ζ(2)Ψ₃(c) − ζ(3)Ψ₂(c) + ζ(4)[Ψ(c) + G(c)] + ζ(2)²G(c) − ζ(5)[Φ(c) + 2F₂(c)] − 2ζ(2)ζ(3)F₂(c),

where:

- **R₅** = [1]M₅ = `hpl_fit_k5_1.json`. It has 310 terms c^eH_w(±c), with 1, 8, 27, 84, 190 terms of weight 2, 3, 4, 5, 6, and e ∈ {−4,…,1}. It has no constant term, and every coefficient has denominator 1, 2 or 4.
- **Ψ₃** = [ζ(2)]M₅ = `hpl_fit_k5_z2.json`. It is the constant ½ plus 60 terms, with 3, 6, 19, 32 terms of weight 1, 2, 3, 4, and e ∈ {−4,…,1}.
- **G** = [ζ(2)²]M₅ = `hpl_fit_k5_z2^2.json` (11 terms and 2 constants). In closed form:
  G = 7/16 − 3/(4c²) − (16 − 19c² + c⁴)/(32c⁴)·L + (7 − 5c²)/(16c³)·A + (c² + 3)(c² − 1)/(32c³)·D.
- **The ladder parts.** They were fitted independently, and they equal the ladder predictions term by term:
  - [ζ(3)]M₅ = −Ψ₂ (`hpl_fit_k5_z3.json`, 22 lines);
  - [ζ(4)]M₅ = Ψ + G (`hpl_fit_k5_z4.json`);
  - [ζ(5)]M₅ = −Φ − 2F₂ = −5/8 + 3/(8c²) − 3(1−c²)²/(8c³)·A (`hpl_fit_k5_z5.json`);
  - [ζ(2)ζ(3)]M₅ = −2F₂ = −1/8 + 3/(8c²) − (1−c²)(3+c²)/(8c³)·A (`hpl_fit_k5_z2*z3.json`).

**With real zeta values,** ζ(4) = π⁴/90 and ζ(2)² = π⁴/36, so the π⁴ part is π⁴(Ψ/90 + 7G/180). Explicitly:

π⁴·[ (65c² − 84)/(2880c²) − (55c⁴ − 165c² + 112)/(5760c⁴)·L − (51c² − 49)/(2880c³)·A + (1 − c²)(3c² − 7)/(1920c³)·D ]

(`s06`, exact through c^200 against the term lists). Here Ψ = ½ + (1/(2c²) − 3/4)L − A/(2c) + (1−c²)/(4c)·D (`s06`, exact).

**Coefficient content (PROVED, `s01`/`s03`).** μ₅,ₙ ∈ ℚ + ℚζ(2) + ℚζ(3) + ℚζ(4) + ℚζ(2)² + ℚζ(5) + ℚζ(2)ζ(3). All seven monomials occur, and no other monomial can, because the weight is at most 5. In Z-words Z_w(2n) with ℚ(n) coefficients:

- the 1-part uses 70 words (depth ≤ 5, indices ±1…±5);
- the ζ(2)-part uses 12 words, the ζ(3)-part 5, and the ζ(4)- and ζ(2)²-parts 2 each;
- [ζ(5)]μ₅,ₙ = −3/((2n−1)(2n+1)(2n+3));
- [ζ(2)ζ(3)]μ₅,ₙ = 2n/((2n−1)(2n+1)(2n+3));
- [ζ(2)²]μ₅,ₙ = −n/((2n−1)(2n+1)(2n+3))·Z₋₁(2n) + (a rational function with denominator 4(n+1)(n+2)(2n−1)(2n+1)²(2n+3)²).

The full formulas are in `s03_prove.out`. The denominators of μ₅,ₙ are only 2n−1, n, 2n+1, n+1, 2n+3 and n+2. These denominators make the range e ∈ {−4,…,1} plausible, but the range was an assumption of the fit. It is not forced by them. What makes it meaningful is this. Harmonic polylogarithms are linearly independent over ℚ(c) (Deneufchâtel, Duchamp, Hoang Ngoc Minh and Solomon, 2011). So the rational function that multiplies each H_w(c) is unique. Hence the prefactors c^e, including the c⁻⁴ ones, are intrinsic to the proved representation. They are not an artefact of the basis, for representations by H_w(±c) with words in 0, ±1 not ending in 0.

## 3. Status of each part

| Part | Terms (lines) | Weights of H_w | Found by | Status |
|---|---|---|---|---|
| R₅ = [1]M₅ | 310 | 2–6 | fit, unique (rank 2186/2186) | **PROVED, all n** |
| Ψ₃ = [ζ2]M₅ | 60 + const ½ | 1–4 | fit, unique (242/242) | **PROVED, all n** |
| G = [ζ2²]M₅ | 11 + 2 const | 1–2 | fit, unique (26/26) | **PROVED, all n** (closed form above) |
| [ζ3]M₅ = −Ψ₂ | 21 + const | 1–3 | fit (80/80) and ladder | **PROVED** (both ways) |
| [ζ4]M₅ = Ψ + G | 11 + 2 const | 1–2 | fit (26/26) and ladder | **PROVED** |
| [ζ5]M₅ = −Φ − 2F₂ | 3 + 2 const | 1 | fit (8/8) and ladder | **PROVED** |
| [ζ2ζ3]M₅ = −2F₂ | 3 + 2 const | 1 | fit (8/8) and ladder | **PROVED** |

"PROVED, all n" means three exact checks for each part separately (`s03_prove.py`). The even identity [c^{2n}](HPL side) − [z]P₅,ₙ = 0 holds word by word as ℚ(n)-combinations of Z_w(2n), for every n ≥ 1. The odd identity [c^{2n+1}] = 0 holds the same way. And the coefficients of c^N for −6 ≤ N ≤ 1 are zero in exact rationals. This is the method and code of G8 (`g8_allN.py`, unchanged), with the left side computed from the definition for symbolic n.

Both sides converge on |c| < 1, so each part is an identity of analytic functions there. As in G8, the theorem is proved for formal coefficients, which implies the real-number formula.

"Unique" has two levels. First, the ranks in the table show it is unique among the terms c^eH_w(±c) with −4 ≤ e ≤ 1 and |w| ≤ 6 − wt(z), plus the constants 1 and c⁻². The referee reproduced these ranks modulo a prime, which implies them over ℚ. Second, by the linear independence of harmonic polylogarithms over ℚ(c) (see §2), each part is also unique as rational-function coefficients of harmonic polylogarithms. It does not mean that fewer terms of another kind cannot do it.

## 4. Ladder predictions vs independent fits (all hold)

There are two independent tests (`s03_prove.out`):

1. **On the μ side, exact for all n.** These are identities of Z-forms computed from the definition, with no fit:
   - [ζ3]P₅ = −[ζ2]P₄;
   - [ζ4]P₅ = [ζ2]P₃ + [ζ2²]P₅;
   - [ζ5]P₅ = −[ζ2]P₂ − 2[ζ2²]P₄;
   - [ζ2ζ3]P₅ = −2[ζ2²]P₄.

   All four are True.
2. **On the function side.** The four term lists fitted without using the ladder are identical, as canonical term lists, to −Ψ₂, Ψ + G, −Φ − 2F₂ and −2F₂. All four are True.

**No ladder prediction failed.**

## 5. Checks

| Check | Label | Result | Command |
|---|---|---|---|
| Symbolic μ₅,ₙ equals the independent routine `g6_common.mu_formal`, n = 1..8, all monomials | exact | True | `python s01_mu_side.py` (~30 s) |
| Fit machinery reproduces G8's R₄ (109 terms) and Ψ₂ (22 lines) exactly | exact | True | `python s00_fit_selftest_k4.py` (~10 s) |
| Fits of all 7 parts (mod 3 primes, CRT plus rational reconstruction, stable from 2 to 3 primes) | fit only | all consistent and unique | `python s02_fit.py 1 z2 'z2^2' z3 z4 z5 'z2*z3'` (~20 s) |
| All-n identities, each part (even, odd, low orders), 4 ladder rungs, 4 fit-vs-ladder comparisons | **PROVED** | all True | `python s03_prove.py` (~1 min) |
| h_w(2n+δ) Z-forms vs the definition: 236 words, δ = −1..5, n = 1..6 (9912 checks) | exact | all equal | `python s04_controls.py` (a few minutes; not timed) |
| Mutants: one coefficient +1, or one exponent −1, in R₅, Ψ₃ and G (6 mutants) | exact | all 6 caught | same |
| Theorem vs definition (Cauchy integral in u, 64 nodes, radius 0.2, 60 dps), at α = π/3, π/2, 2π/3, 1.8, 5π/6, π/4 | NUMERICAL | differences 2.7e-52, 5.5e-53, 6.6e-53, 6.1e-53, 6.7e-53, 3.2e-52 | `python s05_numeric.py A` (~20 s) |
| Theorem at c = 1/√2 vs the stored 400-digit M₅ (`../M5a_next_coefficient/M5_400.txt`) | NUMERICAL | 2.2e-141 (N = 900), 8.1e-170 (N = 1150, at the 170-dps limit) | `python s05_numeric.py B` (~2 s) |
| Closed forms of G, [ζ5], [ζ2ζ3], Ψ and the π⁴ part vs the term lists | exact through c^200 | True | `python s06_small_closed_forms.py` |
| [ζ2^e]P₂ₑ,ₙ for e = 1, 2, 3 (symbolic n) | exact | see §6 | `python s07_top_parts.py` (~4 min) |

All commands run from this folder with `../../../.venv/bin/python`.

The numerical side of A sums each harmonic polylogarithm from its own Taylor series. The definition side does not use the series μ₅,ₙ. Values from the definition:

- M₅(π/3) = 0.008014068707441025871309743475457395331…
- M₅(π/2) = 0.003161648710412257483507878511454280647…
- M₅(2π/3) = 0.0006757075053200854936462481917930191130…
- M₅(1.8) = 0.001762459618682137017208665219691106471…
- M₅(5π/6) = 0.00004063942527594277445117137303918577196…
- M₅(π/4) = 0.01027173693082089109975890223802772617…

## 6. Things found along the way

- **Field rule at k = 5 (PROVED).** The paper already states, in the proof sketch of conj:field, that for k ≤ 7 the denominators are powers of c alone. What is new here is the proved k = 5 representation. Every prefactor is a monomial c^e with −4 ≤ e ≤ 1 (the constants 1 and c⁻² included). So the denominators are powers of c only, consistent with that statement. By the uniqueness remark in §2 these prefactors are canonical. So "the range of e grows by one" is a correct statement about these canonical representations: it was −3…1 at k = 4.
- **Weights.** wt(z) + |w| ≤ 6 holds in every part, with equality reached. R₅'s lower weights 2, 3, 4 have 1, 8, 27 terms, the same counts as R₄'s weights 2, 3, 4. I did not test whether these terms are the same.
- **The π⁴ bracket is no longer elementary at k = 5.** At k = 4, the π⁴ part had only artanh and rational functions. At k = 5, by the ladder, [ζ(4)]M₅ = Ψ + G contains the dilogarithmic D. The real π⁴ coefficient keeps D with coefficient (1−c²)(3c²−7)/(1920c³), which is not zero. So "the π⁴ part is elementary" does not continue. The reason it is not elementary is this. Suppose the bracket were r₁ + r₂·artanh c with r₁, r₂ rational functions of c. Then this would be a nontrivial ℚ(c)-linear relation among H₀,₁, H₀,₋₁, H₋₁,₁, H₁,₋₁ (which occur only in D), H_{±1} and 1. That is impossible, because harmonic polylogarithms are linearly independent over ℚ(c) (Deneufchâtel, Duchamp, Hoang Ngoc Minh and Solomon, 2011) and the coefficient (1−c²)(3c²−7)/(1920c³) is not zero. The same argument shows the bracket is not a polynomial in logarithms either. At the symbol level, the antisymmetric part of S(D) is −(1−c)∧c + (1+c)∧c + 2(1−c)∧(1+c), which is not zero. The referee also verified the explicit π⁴ bracket exactly through c^301.
- **What does continue: the top parts.** These are the parts [z]M_k with wt(z) = k. At k = 5 they are ζ(5) and ζ(2)ζ(3), and they are elementary: rational functions of c and artanh c only. Remark, without proof: by the ladder, every top part is a rational combination of E_e = [ζ(2)^e]M_{2e}. (A one-line weight count from the kernel would settle it. I did not write it out.) We have E₁ = Φ and E₂ = F₂. `s07` gives, for symbolic n the three formulas below. They are exact for all n by computer algebra, as long as that code is right at k = 6. The referee confirmed [ζ2³]μ₆,ₙ independently for n ≤ 80 only. For the other two the referee's exact check covers n ≤ 150.
  - [ζ2]P₂,ₙ = 1/((2n−1)(2n+1));
  - [ζ2²]P₄,ₙ = −n/((2n−1)(2n+1)(2n+3));
  - [ζ2³]P₆,ₙ = n(5n−1)/(3(2n−1)(2n+1)(2n+3)(2n+5)).

  So E₃ is again artanh c plus a Laurent polynomial in c. Its partial fractions only have the denominators 2n+1, 2n+3, 2n+5 and 2n−1, all odd.

  Why the coefficients are rational in n (sketch, not written out): in the ladder's operator form, (∂_a∂_b)^e kills everything in E_{n,2e} except its y^j constants, so no harmonic numbers survive. The pattern of odd denominators only is exact for e ≤ 3 and not proved in general.
- **[ζ(2)²]M₅ carries the Laurent polynomial −3/c³ + 2/c + c** (as the coefficient of D, divided by 32). The same polynomial appears in G2's (3,2) pairing [Y,[X,Y]]⊗X∧Y of M₄. I did not look into whether this is a coincidence.

## 7. What is not done

- **Depth and symbol structure of R₅ were not attempted.** This includes the (3,3)/(4,2) cobracket of its weight-6 part and the point-level questions.
- **The term lists are data.** They are unique in the stated basis, but I did not try to shorten them, for example into classical polylogarithms or with fewer pairs.
- **Control.** The all-n proof reuses G8's code unchanged; the new pieces are the k = 5 driver `s03_prove.py`, the fit `m5fit.py` and the numerics `s05_numeric.py`. G8's own k = 4 output is reproduced by `s00` as a control.

## 7a. Referee

This is an independent AI re-check with its own code (`new_results/M5AR_referee_M5_every_angle/REFEREE.md`). It is not peer review. Verdict: confirmed.

- The referee's own evaluation of M₅ from the definition agrees with the closed form to 1.3e-98 at eight angles.
- Every part agrees exactly for all n ≤ 150.
- 12 of 12 planted errors were caught.
- The referee found no error in any formula, coefficient, count or number. It asked for the rewordings that are now in this file.

What the referee did not check:

- The written proofs of the general-k lemmas behind `g8_allN.mu_side` (Lemmas g6:lem:En and g6:lem:A, Definition g6:def:formal). The k = 5 μ-side for every n rests on them. The independent agreement covers n ≤ 150, and n ≤ 80 at k = 6.
- `g8_zsum.py` line by line. Its k = 5 use is covered by the exact comparison for n ≤ 150 and by the 9912-check selftest (s04, not rerun).
- It did not rerun s00, s02, s04, s05, s06 and s07. It re-derived their substance with its own code, except the k = 6 computation (n ≤ 80 only).
- Depth, symbols or cobrackets of R₅.
- Theorem thm:zladder in general. Only the four k = 5 rungs and [ζ3]M_{k+1} = −[ζ2]M_k for k ≤ 4.

## 8. Files

- `result.md` — this report.
- `s01_mu_side.py` / `.out` — symbolic μ₅,ₙ, saved to `mu5_zsum.pkl`.
- `m5fit.py`, `s00_fit_selftest_k4.py` / `.out`, `s02_fit.py`, `s02_fit_small.out`, `s02_fit_R5.out` — the fit (not a proof).
- `hpl_fit_k5_{1,z2,z2^2,z3,z4,z5,z2*z3}.json` — the term lists. Quote the names in a shell.
- `s03_prove.py` / `.out` — the all-n proof, the ladder rungs, and the fit-vs-ladder comparison.
- `s04_controls.py` / `.out` — the machinery self-test and the mutants.
- `s05_numeric.py`, `s05_numeric_A.out`, `s05_numeric_B.out` — numerical checks.
- `s06_small_closed_forms.py` / `.out` — closed forms of G, [ζ5], [ζ2ζ3], Ψ, and the π⁴ part.
- `s07_top_parts.py` / `.out` — the top parts E₁, E₂, E₃.
- Copied, unchanged, from `paper/ieee/proofs/scripts/`: `g8_allN.py`, `g8_zsum.py`, `g6_common.py`, `g8_terms_k3_partz2.json`, `g8_terms_k4_partz2.json`, `g8_terms_k4_part1.json`.
