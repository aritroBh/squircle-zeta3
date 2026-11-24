# CHN1: Irreducible MZVs in the volume of the chain body B(P4)

**Work directory:** `<project>/new_results/CHN1_chain_bodies/` (all code + checkpoints there; step files `step_01.py`–`step_87.py`).

**Object.** B = {x ∈ ℝ⁴ : |x₁|^p₁+|x₂|^p₂ ≤ 1, |x₂|^p₂+|x₃|^p₃ ≤ 1, |x₃|^p₃+|x₄|^p₄ ≤ 1}, uᵢ = 1/pᵢ = aᵢu. Vol/16 = G·₃F₂(−u₄,u₃,1+u₁+u₃; 1+u₃,1+u₁+u₂+u₃; 1), G = Γ(1+u₂)Γ(1+u₁+u₃)/Γ(1+u₁+u₂+u₃). Everything below was computed by an exact engine (Fraction arithmetic; ₃F₂ expanded into stuffle harmonic sums, summed to MZV words, reduced by the finite double-shuffle + Hoffman regularized relations) and cross-checked numerically.

---

## (i) The ζ(3,5) form P8 and the quartic Q4

c₈(a) = [u⁸]Vol/16 lies in span{ζ(8), ζ(3)ζ(5), ζ(2)ζ(3)², ζ(3,5)}; write its ζ(3,5)-coefficient as P8(a). **Exact result:**

```
P8(a) = a1 a2 a3 a4 (a1+a2+a3+a4) C3(a)/10
Q4(a) = P8/(a1a2a3a4) = (a1+a2+a3+a4) C3(a)/10
```

with (x,y,z,w) = (a₁,a₂,a₃,a₄):

```
C3 = 4(x³+w³) − 4(y³+z³) + 3wx(x+w) − 3yz(y+z)
     + 6(x²+w²)(y+z) − 6(y²+z²)(x+w) + 6yz(x+w) − 6wx(y+z)
   = 2(p−q)(2p+q)(p+2q) − 9·a1a4·(p+2q) + 9·a2a3·(2p+q),    p = a1+a4, q = a2+a3.
```

**Factorisation over ℚ:** (a₁+a₂+a₃+a₄) × C3, with C3 irreducible over ℚ; C3 is a **smooth cubic surface** in P³ (Gröbner basis of (C3, ∇C3) has only the trivial zero ⇒ smooth ⇒ absolutely irreducible). It contains the projective lines Π1, Π2 (below).

**Vanishing tests for Q4 (exact polynomial identities):**

| locus | Q4 ≡ 0? | restriction of C3 |
|---|---|---|
| Π1 = {a₁=a₂, a₃=a₄} | **YES** | 0 |
| Π2 = {a₁=a₃, a₂=a₄} | **YES** | 0 |
| diagonal (t,t,t,t) | point zero | C3(1,1,1,1)=0 |
| plane a₁+a₄ = a₂+a₃ | no | 27(a₁−a₂)(a₁−a₃)(a₂+a₃) — on this plane the positive zero set is exactly Π1 ∪ Π2 |
| palindromic plane (x,y,y,x) | no | 2(s−t)(7s²+13st+7t²) (positive zero only s=t) |
| locus a₁=a₂=a₃ | no | −(s−t)(4s²+19st+4t²) (positive zero only s=t) |

**Gradient data (reversal + Euler force ∝(1,−1,−1,1)):** ∇C3(1,1,1,1) = 27(1,−1,−1,1); ∇Q4(1,1,1,1) = (54/5)(1,−1,−1,1); ∇P8(1,1,1,1) = (54/5)(1,−1,−1,1). Neither factor of Q4 vanishes at the diagonal and the gradient is nonzero ⇒ **the diagonal is a smooth point of the surface {Q4 = 0}** (the tangent plane there is a₁−a₂−a₃+a₄ = 0, which meets the surface in the three lines Π1, Π2, and {a₂=a₃=0, a₁+a₄=0}).

**ζ(3,5) coefficients (engine, exact):**

| direction a | P8(a) | Q4(a) |
|---|---|---|
| (1,1,1,1) | 0 | 0 |
| (1,1,1,2) | 58 | 29 |
| (2,1,1,1) (reversal) | 58 | 29 |
| (1,2,2,1) | −1464/5 | −366/5 |
| (3,1,2,5) | 36696 | 30616/5 |
| (5,2,1,3) (reversal) | 36696 | 30616/5 |
| **(3,1,1,1) blind** | **1746/5** | **582/5** |

P8 = a₁a₂a₃a₄·Q4 is reversal-invariant; the trivial divisibility by a₃a₄ (and the a₁a₂ factor via the Gauss boundaries, control E) is confirmed; Q4 was also recovered independently by interpolating the 19 reversal-invariant quartic monomials from 19 engine points, validated at 21 held-out points (exact).

---

## (ii) Weight-10 and weight-11 irreducible coefficient forms, factored

All forms divisible by a₁a₂a₃a₄. Σa = a₁+a₂+a₃+a₄; C3, Q5 as given (Q5 = the 56-term quintic printed in step_86, irreducible over ℚ).

| basis element | coefficient form (exact) | zero set on positive directions |
|---|---|---|
| ζ(2)ζ(3,5) (w=10) | −a₁a₂a₃a₄(a₂+a₄)(a₁+a₃)Σa·C3/10 | = {C3=0} (other factors >0) |
| ζ(3,7) (w=10) | a₁a₂a₃a₄·Σa·Q5/224 | ⊇ Π1 ∪ Π2 ∪ {diagonal} |
| ζ(3)ζ(3,5) (w=11) | a₁a₂a₃a₄(a₂+a₄)(a₁+a₃)(Σa)²·C3/10 | = {C3=0} |
| ζ(3,3,5) (w=11) | −a₁a₂a₃a₄·Σa·mid·C3/10, mid = w²x+wx²+2wxy+2wxz+2wyz+2xyz+y²z+yz² (irreducible, positive on positive directions) | = {C3=0} |

So all ζ(3,5)-type forms vanish exactly on {C3 = 0}; the ζ(3,7) form has its own quintic Q5, which **vanishes identically on Π1, on Π2, and at the diagonal** (Q5(1,1,1,1)=0, ∇Q5(1,1,1,1) = 2988(1,−1,−1,1)), but **not** on {C3=0} generally: at the positive point (1/2, 2, 2, w≈2.3995) ∈ {C3=0} (off Π1∪Π2), the ζ(3,5)-forms are 0 while the ζ(3,7) form ≈ 1.62e3 ≠ 0. **On the zero set of Q4 outside Π1∪Π2, only ζ(3,5) (and its products) vanishes; ζ(3,7) does not.**

**Values of the forms:**

| a | ζ(2)ζ(3,5) | ζ(3,7) | ζ(3)ζ(3,5) | ζ(3,3,5) |
|---|---|---|---|---|
| (1,1,1,1) | 0 | 0 | 0 | 0 |
| (1,1,1,2) | −348 | 29515/56 | 1740 | −1276 |
| (1,2,2,1) | 13176/5 | −45849/14 | −79056/5 | 61488/5 |
| (3,1,2,5) | −1100880 | 21250515/14 | 12109680 | −9100608 |
| (3,1,1,1) blind | −13968/5 | 142173/28 | 83808/5 | −59364/5 |

(At weights 9 and 11 no other irreducible directions occur: the quotient has dimensions 5 and 9 with only ζ(3,3,5), ζ(3)ζ(3,5) involving irreducibles.)

---

## (iii) The exponent rule

**Theorem (EXACT, by computer algebra over ℚ):** through weight 11, no irreducible MZV and no product containing one appears in any coefficient c_k(a) (k ≤ 11) **if and only if** a ∈ Π1 = {a₁=a₂, a₃=a₄} or a ∈ Π2 = {a₁=a₃, a₂=a₄} (among positive directions; the two planes meet at the diagonal a₁=a₂=a₃=a₄).

- *If:* all five forms above vanish identically on Π1 and on Π2 (polynomial identities: C3 contains both lines; Q5 vanishes on both).
- *Only if:* the common zero set of the forms is {C3=0}∩{Q5=0}. The resultant of C3 and Q5 (a₃=1) in a₄ factors completely into 14 lines: 729000000·x·y·(x−1)(x+1)(x+2)(x−y)(x+y)(x+2y)(2x+1)(2x+y)(y+1)(x−y+1)(x+y−1)(2x+y+1); solving for a₄ over each shows **all components except Π1 (a₁=a₂, a₄=a₃) and Π2 (a₁=a₃, a₄=a₂) have a nonpositive coordinate**. (The a₁=0, a₂=0 boundaries are the Gauss controls, not part of the statement.)

**Through weight 12 the rule still holds**, but by a new mechanism (see §(v) below): on Π1, Π2, and the diagonal, c₁₂ has *nonzero* ζ(3,9) and ζ(5,7) coefficients in exact ratio −14/9, so the combination (14/9)ζ(3,9)+ζ(5,7)-part collapses into products via the weight-12 relation; the coefficients of the exceptional generator ζ(2,1,1,8) and of ζ(2)²ζ(3,5), ζ(2)ζ(3,7) vanish identically on both planes (proven at 15 exact rational points per plane; the diagonal lies on both). So c₁₂|planes ∈ span(products) as well.

**All-orders identity: NOT found (precise negatives).**
1. On each of Π1, Π2, diagonal, log Vol/16 = (single-zeta series in the locus parameters) + (a product series **identical on all three loci**, in the DS quotient, through weight 11). If the volume were (Γ-product)×(universal factor), the single-zeta parts would be power sums Σ±ℓᵢ(x,y)^m of linear forms. An exact solver over forms αs+βt, |α|,|β|≤3, reproduces Q_m for m ≤ 8 only with huge non-integral multiplicities and **fails from m = 9** (e.g. Π2: Q₂ = −6st, Q₃ = −15st(s+t), Q₅ = −85st(s+t)(s²+st+t²); diagonal: Q_m/s^m = −6, −30, −76, −510, −954, −8190, −37448/3, …, with odd m fitting 2−2^{2m−1}, i.e. a fractional power Γ(1+4u)^{-1/2}-type obstruction). So **no single Γ-product times a universal factor** (limited-search evidence).
2. Two-term ansatz Vol|Π1 = R₁(s,t)G(s,t)+R₂(s,t)G(t,s) with G the natural Γ-product Γ(1+s)Γ(1+s+t)/Γ(1+2s+t): **fails already at weight 4** (exact; R₁ = s²/(s²−t²), R₂ = −t²/(s²−t²) forced by weights 2,3, then weight 4 gives 5(s+t)²(s²+t²)/4 ≠ 0).
3. **Thomae S₅-orbit (all 120 elements computed exactly):** on Π1 and on Π2 the ₃F₂ is Thomae-equivalent to a **Saalschützian (σ=1) ₃F₂**, ₃F₂(s, t, 1+2s+2t; 1+2s+t, 1+s+2t; 1) in the plane's two parameters (note 1+2s+2t = c, denominators c−s, c−t) — but **no orbit element is Gauss-collapsible, Dixon-, Watson-, or Whipple-summable**, and **the diagonal's orbit contains no Saalschützian element at all**. So no classical summation explains the vanishing; the Saalschützian form is a reformulation, not an evaluation.
4. Step 9(c) (Schlotterer–Stieberger M3/M5 commutation on the planes): **not done** (no literature check was performed).

---

## (iv) Weight 12 (optional step 10)

- DS+Hoffman table at weight 12: 1024 admissible words, **nullity 12** (conjectural value); duality holds for all 1024 words; sum theorem holds.
- **Generators chosen:** {ζ(12), ζ(3)ζ(9), ζ(2)ζ(3)ζ(7), ζ(2)ζ(5)², ζ(2)²ζ(3)ζ(5), ζ(2)³ζ(3)², ζ(3)⁴, ζ(2)²ζ(3,5), ζ(2)ζ(3,7), ζ(3,9), ζ(5,7), ζ(2,1,1,8)}. Products span dimension 10; the 12th direction is the exceptional word ζ(2,1,1,8) (dual of ζ(2,1⁶,4)); ζ(3,3,6), ζ(3,4,5), ζ(2,3,7), ζ(3,3,3,3) all lie in the products+ζ(3,9)+ζ(5,7) span.
- **The Broadhurst–Kreimer-type relation (exact in the quotient, and validated numerically to 6.2e-61):**
  ζ(5,7) = (128195/12438)ζ(12) − (14/9)ζ(3)ζ(9) − (25/3)ζ(5)ζ(7) + (14/9)ζ(3,9), i.e. ζ(5,7) ≡ (14/9)ζ(3,9) mod products.
- Coefficients at generic points: (1,1,1,2): ζ(3,9) = 418268281/900, ζ(5,7) = −44869481/150, w* = −9280/3, ζ(2)²ζ(3,5) = −8236/5, ζ(2)ζ(3,7) = −4695/4. Blind (3,1,1,1): 16766672/5, −21333369/10, −24600, −540792/25, −173646/7. Diagonal: ζ(3,9) = 41328, ζ(5,7) = −26568 (exact ratio −14/9), w* = 0. Diagonal c₁₂ validated by an independent Cauchy/FFT computation: |diff| = 4.8e-71.

---

## (v) Labels

- **EXACT (computer algebra over ℚ, not a formal proof):** all coefficient forms and factorisations ((i)–(iv)), the vanishing identities on Π1/Π2/diagonal, the only-if component analysis, the weight-12 relation, and the "no irreducibles through weight 12 on Π1∪Π2" statement. These are identities of real numbers provided the engine is correct (supported by controls A–E, the pilot-vector reproductions, the pilot-file agreement through c₁₁, and two independent Cauchy checks; **not certified**) and every reduction-table entry came from proved relations (finite double shuffle + Hoffman regularization). **No PSLQ was used anywhere.**
- **CONJECTURAL:** ℚ-linear independence of the basis elements (only the dimension upper bounds are theorems; our nullities match the conjectured dimensions 1,1,1,2,2,3,4,5,7,9,12). This is what makes "ζ(3,5) is genuinely needed" true for real numbers. The direction "coefficient = 0 ⇒ only single-zeta products" is unconditional given the table.
- **KNOWN:** the MZV mechanism (₃F₂(1) and M_{0,5} expansions: Terasoma 2002; Schlotterer–Stieberger arXiv:1205.1516; Broedel–Schlotterer–Stieberger–Terasoma arXiv:1304.7304). **No literature check for the novelty of P8/Q4 was done.**
- **NUMERICAL:** the blind check (c₈ at (3,1,1,1)) and the diagonal c₁₂ Cauchy check.

## (vi) Failed controls and unexpected nonzeros

- **No control failed.** A–E all passed: exact diagonal c₂..c₆, 3c₈; the four pilot c₈ vectors (3c₈, 32c₈, 80c₈, 64c₈); reversal invariance of all weight ≤ 8 coefficient vectors at 6 direction pairs; Gauss Γ-products at a₁=0 (Γ(1+u₃)Γ(1+u₂+u₄)/Γ(1+u₂+u₃+u₄)) and a₂=0 (Γ(1+u₃)Γ(1+u₄)/Γ(1+u₃+u₄)) matching whole vectors at 4+4 directions; pilot diagonal 5c₁₀ vector.
- **Unexpected nonzero (resolved):** at weight 12, ζ(3,9) and ζ(5,7) DO have nonzero coefficients on Π1, Π2, and the diagonal — but in exact ratio −14/9, so they cancel modulo the weight-12 relation back into products. First flagged here; this is the precise sense in which the "no irreducibles" phenomenon survives weight 12.
- The pilot file `chain4_generic_3_1_2_5.out` line 1 (stale nested-quadrature artefact) was ignored per instructions.
- **Reported negatives:** singles are not small power sums (no Γ×universal identity); two-term Γ ansatz fails at weight 4; no classically summable Thomae image; the diagonal's Thomae orbit is not even Saalschützian; step 9(c) not done.

## Evidence table (residuals and digit counts)

| check | method | residual | digits |
|---|---|---|---|
| G·₃F₂ vs x-coordinate integral, 3 rational points | mpmath quad vs hyp3f2, dps 35 | 7.5e-37 / 0 / 1.9e-32 | ~30 |
| Vol(uᵢ=1/2)·16 = 4+π²/2 | hyp3f2 dps 35 | 0 | 35 |
| anchors ζ(3,5), ζ(5,3) | mzv_vals.dz | matches to given 20 digits | 20+ |
| 6 depth-2 weight-8 words vs dz(3,5)-based reduction | dps 80 | ≤ 1.8e-80 | ~80 |
| 5 all-≥2 words (5,3),(2,3,3),(3,2,3),(4,2,2),(2,2,2,2) Richardson vs DS | nested sums + Richardson, dps 80 | 8.8e-49 … 1.0e-48 | ≥ 48 |
| exact c₂..c₈ vs 6 pilot files | dps 120 evaluation | diag 4.6e-51..3.0e-49 (tol 1e-45); generic 1.2e-70..4.6e-65 (tol 1e-60) | ≥ 60 |
| exact c₉,c₁₀,c₁₁ vs pilot files | dps 90 | c₉,c₁₀: ≤ 4.9e-66; c₁₁: within pilot aliasing (~1e-41) | as pilot allows |
| **blind check c₈ at (3,1,1,1)** | Cauchy/FFT N=128, r=R/6=1/30, dps 85, u₃∫₂F₁ representation | **4.23e-71** (imag −2.5e-76, c₀ = 1) | 71 |
| diagonal c₁₂ | Cauchy/FFT N=128, r=1/18, dps 85 | 4.84e-71 (imag 2.2e-72) | 71 |
| weight-12 relation ζ(5,7)=… | dz numerics dps 60 | 6.17e-61 | 61 |
| Q4 interpolation (19 unknowns) | exact solve + 21 held-out points | exact | exact |
| symbolic engine vs single-direction | 6 random directions, all W ≤ 11 | exact | exact |

(Note on the blind-check settings: the task's suggested N=64 at r=R/6 gives FFT aliasing ≈ (1/6)⁶⁴·R⁻⁸ ≈ 4e-48 for R≤1/2-class directions, insufficient for the 1e-60 target; I used the pilot-proven N=128 with r=R/6, dps 85, aliasing ~1e-94. The ₃F₂ was evaluated via the analytically continued identity ₃F₂ = 1 + u₃∫₀¹ t^{u₃−1}(₂F₁(t)−1)dt, verified against mpmath hyp3f2 to 7.9e-89, ~1s/point vs ~72s. Prediction was written to `blind_prediction.txt` before computing.)

## Files (all in the work directory) and reruns

- `mzv_exact.py` — words, stuffle/shuffle, DS+Hoffman relations, pure-ℚ and flint rref, reduction tables, duality. Test: `python mzv_exact.py`
- `mzv_num.py` — Richardson nested-sum MZV evaluator.
- `engine.py` — single-direction exact engine (c_W(a) → DS basis). Controls: step_18–step_21 (run via `python -c "import engine; print(engine.reduced_coeffs((1,1,1,1),8))"`).
- `symengine.py` — symbolic engine in u₁..u₄, total degree ≤ 11 → `forms.pkl` (all coefficient forms).
- `mzv_vals.py` (copied from pilot) — dz(a,b).
- `blind_check.py` — blind c₈ job → `blind_points.jsonl`, `blind_result.json`; prediction `blind_prediction.txt`.
- `w12_planes.py` — 15-point/plane weight-12 vanishing proof → `w12_planes.json`.
- `diag_c12.py` — diagonal c₁₂ Cauchy check → `diag_c12_result.json` (+ `diag_c12_prediction.txt`).
- `locus_logV.pkl` — cached log Vol coefficients on the three loci.
- `progress.log` — one line per finished step; `transcript.json` — full harness log.
- Rerun the main computations: the harness step files `step_09.py`–`step_87.py` (each is self-contained; earlier ones build `mzv_exact.py`/`engine.py`/`symengine.py`, later ones the analysis). E.g. `python step_24.py` (symbolic engine), `python step_36.py` (Q4 interpolation), `python step_57.py` (W=12 table), `nohup nice -n 10 python blind_check.py &`.

**Bottom line.** The exponent vectors that admit irreducible MZVs into Vol B(P4) are exactly (positive directions) the complement of Π1 = {a₁=a₂, a₃=a₄} ∪ Π2 = {a₁=a₃, a₂=a₄}, through weight 11 — and still through weight 12, the planes' escape at weight 12 being rescued only by the weight-12 relation ζ(5,7) ≡ (14/9)ζ(3,9) mod products. The weight-8 irreducible enters via P8 = a₁a₂a₃a₄(a₁+a₂+a₃+a₄)·C3/10 with C3 a smooth cubic surface containing the two plane-lines and passing smoothly through the diagonal; ζ(3,7), ζ(3,3,5), and the ζ(3,5)-products enter through the factored forms tabulated above. No all-orders Γ-product or classical-summation identity explaining the planes was found; the Saalschützian Thomae image on Π1/Π2 is the remaining lead.