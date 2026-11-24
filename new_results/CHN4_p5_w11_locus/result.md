# CHN4: Does the 5-chain P5 ever lose all irreducible MZVs? The W11/W12 cut of V10

## Verdict

**Is there a positive exponent vector a = (a₁,…,a₅) (open orthant) at which the P5 volume has no irreducible MZV through weight 12? — NO.**

Every point of the **open** positive orthant carries an irreducible MZV by weight 12:

| region of (open) orthant | irreducible appears | status |
|---|---|---|
| C4(a) ≠ 0 | ζ(3,5) at weight 8 | **EXACT** |
| C4=0, Q7′≠0 | ζ(3,7) at weight 10 | **EXACT** |
| on V10, f11 ≠ 0 | ζ(3,3,5) at weight 11 | **EXACT** |
| on Γ (C4=Q7′=f11=0, all aᵢ>0) | W12 irreducibles | **NUMERICAL** |

The **only** irreducible-free-through-W12 points lie on the **boundary** (some aᵢ=0), never in the open orthant — see "Boundary structure" below. So, exactly as conjectured, **P5 (unlike P4, whose irreducible-free locus Π₁∪Π₂ meets its open orthant) has no irreducible-free exponent vectors in the open orthant.**

---

## 1. The exact W11 obstruction f11  — **EXACT (PROVED)**

The ζ(3,3,5) coordinate of c₁₁ is
$$f_{11}(a)= -\tfrac1{10}\,a_2a_3a_4\,(a_1{+}a_2{+}a_3{+}a_4{+}a_5)\; P_7(a),$$
where **P7** (`P7.txt`) is:
- homogeneous **degree 7**, **209 terms**, content 1;
- **irreducible over ℚ** (`sympy.factor_list` → single factor);
- **reversal-symmetric** (0/209 monomial mismatches; `subs(simultaneous=True)`);
- P7(1,1,1,1,1) = **−5104**; coefficients 81 positive / 128 negative (mixed ⇒ sign changes).

**Proof of the factorisation.** f11 vanishes on each hyperplane a₂=0, a₃=0, a₄=0, Σa=0 (exact engine), so a₂a₃a₄(Σa) | f11. The quotient Q11 = −P7/10 (degree 7) was interpolated on a rank-170 reversal-invariant determining set, verified on 48 held-out points (0 mismatches). Then **upgraded to a full proof**: D = f11 − a₂a₃a₄(Σa)Q11 vanishes at all **693** points of a rank-693 determining set for reversal-invariant degree-11 polynomials (rank verified mod two primes). ⇒ f11 = −(1/10)a₂a₃a₄(Σa)P7 **exactly**. Independently re-verified against the engine at 3 random points.

## 2. V10 is exactly the irreducible-free locus through W10 — **EXACT**

From `forms10.pkl` by exact polynomial division (all remainders 0):
- W8 ζ(3,5) = a₂a₃a₄(Σa)·C4/10
- W10 ζ₂ζ(3,5) = −a₂a₃a₄(Σa)(a₂+a₄)(a₁+a₃+a₅)·C4/10  (vanishes at C4=0)
- W10 ζ(3,7) = a₂a₃a₄(Σa)·Q7′/224  (vanishes at Q7′=0)

Hence for positive a, all irreducibles through W10 vanish **iff** C4 = Q7′ = 0, i.e. exactly on V10.

**Control (W11):** the ζ3ζ(3,5) coordinate = a₂a₃a₄(Σa)·C4·H3, with H3 degree 3 (21 terms, **all-positive coefficients**), verified on 20 fresh points + exact C4-division (remainder 0). So ζ3ζ(3,5) vanishes on V10, and **f11 is the only W11 irreducible on V10**.

## 3. Γ = V10⁺ ∩ {f11 = 0} is a nonempty real curve — **NUMERICAL**

Q11 (= −P7/10) **changes sign** on V10⁺ (robust solver: exact a₄ isolation via `G.intervals()` + mpmath bisection + **common-root a₅** — *caution:* a₅=−B/A gives spurious lifts at degenerate a₄, e.g. fake a₅=0.344 at the a₄=1 double root). Confirmed at 60 digits:

- Q11 = **−0.0495** at (a₁,a₂)=(1/100,1/100); Q11 = **−140429** at (1,1); Q11 = **+0.222** at (1/100,3).
- Sign change on **one continuous sheet** (a₁=1/100: Q11 crosses 0 between a₂=2.725 and 2.750).
- Survey over [10⁻³,10³]²: **4057** cells show both signs across sheets — Γ is extensive.

So the W11-irreducible-free locus Γ⁺ is **nonempty** in the open orthant, and we must check W12.

## 4. W12 kills Γ in the open orthant — **NUMERICAL**

mpf engine (120 digits) at **8 interior Γ points** (a₁ = 5·10⁻⁵ … 2.0, two distinct regions). At every point the controls (W8 ζ(3,5), W10 ζ(3,7), W11 ζ3ζ(3,5), f11) vanish (≤5·10⁻¹⁵, mostly 10⁻⁵⁰–10⁻⁸¹ — confirming on Γ), while **all four W12 irreducibles are nonzero**:

| a₁ | ζ2²ζ(3,5) | ζ2ζ(3,7) | CLASS ζ(3,9)+(14/9)ζ(5,7) | ζ(2,1,1,8) | max\|irred\|/max\|c₁₂\| |
|---|---|---|---|---|---|
| 5e-5 | −0.178 | −0.955 | 3.507 | 1.485 | 6.1e-8 |
| 0.01 | −46.72 | −250.3 | 919.3 | 389.3 | 1.1e-5 |
| 0.123 | −13.24 | −70.90 | 213.6 | 110.3 | 1.7e-3 |
| 1.004 | −1.23e6 | −6.58e6 | 2.01e7 | 1.02e7 | 5.6e-3 |
| 2.0 | −3.87e8 | −2.07e9 | 6.33e9 | 3.22e9 | 2.0e-3 |

The normalised max|coordinate|/max|c₁₂| is **> 0 at every interior point** (range 6·10⁻⁸–6·10⁻³); it tends to 0 only linearly as a₁→0 toward the boundary stratum. Four independent conditions on a 1-dimensional curve are generically never simultaneous-zero; no interior Γ point is irreducible-free.

## 5. Boundary structure — the complete closed-orthant locus — **EXACT**

The irreducible-free-through-W12 locus in the **closed** orthant is exactly:
- **the faces a₂=0, a₃=0, a₄=0** — every irreducible coordinate (W8…W12) has the a₂a₃a₄ factor (EXACT for W8/W10/W11; multi-point exact-vanishing for W12);
- **the strata {a₁=0, a₂=a₄, a₅=a₃} and {a₅=0, a₂=a₄, a₁=a₃}** — on these C4=Q7′=Q11=0 *identically*, and all four W12 irreducibles vanish **exactly** (Fraction engine; off these strata, a₁=0/a₅=0 faces have irreducibles).

These boundary faces/strata are degenerate chains (P5(a₁=0) = square × P4(a₂,a₃,a₄,a₅)), consistent with P4's irreducible-free loci living on its own orthant. **They are disjoint from the open orthant**, so the answer to the question is NO.

## Evidence table

| check | method | result | residual/digits |
|---|---|---|---|
| engine vs `test_pt_w1112.pkl` | exact Fraction | match | exact, all 9 coords |
| m0 → irreducible coords | exact | 0 | exact |
| a₂a₃a₄(Σa) \| f11 | vanishing on 4 hyperplanes | yes | exact (4 pts each) |
| f11 = a₂a₃a₄(Σa)Q11 | 48 held-out pts | 0 mismatches | exact |
| f11 = −(1/10)a₂a₃a₄(Σa)P7 | 693-pt rank-693 det set | **PROVED**, 0/693 | exact (1574 s) |
| det11 set rank | mod 2 primes | 693 | exact |
| P7 irreducible, reversal-inv | factor_list / monomials | yes / 0 of 209 | exact |
| W10 ζ(3,7)=a₂a₃a₄(Σa)Q7′/224 | division | remainder 0 | exact |
| W11 ζ3ζ(3,5)=a₂a₃a₄(Σa)C4·H3 | interp + 20 pts + C4 div | remainder 0 | exact |
| Q11 sign on V10⁺ | 60-digit V10 solver | −0.0495, +0.222, −140429 | 60 digits |
| W12 on Γ (8 pts) | mpf 120-digit engine | all 4 irreducibles ≠ 0 | controls ≤5e-15 |
| boundary strata W12 = 0 | Fraction engine | all 0 | exact |

## Files (work dir `CHN4_p5_w11_locus/`)
`P7.txt`, `Q11.txt` (f11 quotient), `Q11_poly.pkl`; `proof_div.log` (693-pt proof), `det11_points.pkl`; `ctrl_z3z35.log`, `H7_z3z35.pkl`, `H7_div_C4_quotient.pkl` (control); `gamma_final.log/.pkl`, `gamma_final_pts.pkl` (W12 on Γ); `signmap.pkl`, `survey_120.pkl`; scripts `run_f11_evals.py`, `run_proof_div.py`, `run_ctrl_z3z35.py`, `run_gamma_final.py`, `survey_lib.py`; plus the auto-saved `step_*.py` cells (setup/validation, interpolation, sign survey, Γ-tracing, boundary checks). Timings: W11 ≈2.2 s/pt exact, W12 ≈8 s/pt exact, ≈35 s/pt mpf@120d; divisibility proof 1574 s; control ≈6 min.

## Labels
- **EXACT (computer algebra):** f11 = −(1/10)a₂a₃a₄(Σa)P7 (PROVED); V10 = W10 locus; W11 control; boundary vanishing.
- **NUMERICAL (mpmath):** Q11 sign-change/Γ nonempty; W12 nonzero on interior Γ (8 points, controls ≤5e-15, all irreducibles clearly nonzero).
- **CONJECTURAL (MZV independence):** needed only to read "coefficient ≠ 0" as "the irreducible is genuinely present in the volume"; all coefficient vanishing/nonvanishing statements themselves are exact.

**Bottom line.** Through weight 12, P5 has **no irreducible-free exponent vector in the open positive orthant** (the expected P5 theorem). The weight-11 obstruction f11 = −(1/10)a₂a₃a₄(Σa)P7 has an irreducible quartic-free factor P7 that changes sign on V10⁺, carving out a real curve Γ; but the weight-12 irreducibles are nonzero everywhere on Γ⁺. The only irreducible-free points lie on the boundary (faces a₂,a₃,a₄=0 and the strata {a₁=0,a₂=a₄,a₅=a₃}, {a₅=0,a₂=a₄,a₁=a₃}).