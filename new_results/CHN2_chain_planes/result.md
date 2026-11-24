# CHN2 Part B/C — results

**Work directory:** `<project>/new_results/CHN2_chain_planes/`. All steps auto-saved as `step_01.py`–`step_75.py` by the harness; helper modules written by me: `engine5.py`, `kimilib.py`; engine copies from CHN1 (read-only sources, copied): `mzv_exact.py`, `engine.py`, `symengine.py`, `mzv_num.py`, `mzv_vals.py`.

---

## Part B — cross-check of the planner's closed form on Π1 ∪ Π2: FULL AGREEMENT

Expanded `(B+G0)/2` exactly in u on Π1 = (x,x,y,y) and Π2 = (x,y,x,y) using log Γ(1+x) = −γx + Σ(−1)^m ζ(m)x^m/m (Fraction arithmetic; the γ-terms cancel since all Γ-arguments are balanced), reduced each coefficient to the CHN1 double-shuffle basis, and compared with CHN1's engine:

| weight | method | result |
|---|---|---|
| 2–11 | exact polynomial identity vs `forms.pkl` restricted to each plane, all basis coordinates, both planes | **MATCH (EXACT)** |
| 12 | vs `engine.reduced_coeffs` at 20 rational points per plane, all 12 basis coordinates (degree-12 homogeneous polys in (x,y) have 13 coefficients; 20 ratios decide identity) | **MATCH (EXACT)** |

Weight-12 detail: on both planes the coordinates of ζ(2)²ζ(3,5), ζ(2)ζ(3,7), ζ(2,1,1,8) vanish identically; ζ(3,9) and ζ(5,7) appear with coefficient ratio exactly z₃₉ = −(14/9)z₅₇ (the CHN1 "weight-12 rescue"), matching the engine at all 20 points per plane. **The planner's theorem Vol/16 = (B+G0)/2 on Π1 ∪ Π2 is confirmed through weight 12 with zero mismatches.** (One self-correction en route: my first Π2 substitution of G0 was wrong — u1+u3 = 2xu on Π2, not (x+y)u; caught because Π1 matched and Π2 didn't; the planner theorem was never in doubt.)

File: `partB_planecoords.pkl` (symbolic plane coordinates, W≤12). Rerun: steps 13–18 (see `step_13.py`–`step_18.py`; core machinery in `step_13.py`+`step_17.py`).

---

## Part C — the 5-chain P5

### C0. Volume formulas (PROVED derivations; NUMERICAL checks all passed)

With u_i = 1/p_i = a_i u, t_i = |x_i|^{p_i} (density u_i t_i^{u_i−1}):

1. **t-integral:** Vol/32 = u2u3u4 ∭_{t2+t3≤1, t3+t4≤1} t2^{u2−1} t3^{u3−1} t4^{u4−1} (1−t2)^{u1} (1−t4)^{u5} dt2 dt3 dt4 (integrate out t1, t5).
2. **Double series:** Vol/32 = Σ_{n,m≥0} [(−u1)_n/n!·u2/(u2+n)] [(−u5)_m/m!·u4/(u4+m)] · Γ(1+u3)Γ(1+u2+u4+n+m)/Γ(1+u2+u3+u4+n+m) — note the denominator gamma contains u2+u3+u4, **no u1** (I initially mis-indexed this; the u=(1,1,1,1,1) anchor caught it).
3. **Single series + reversal (the engine form):** Vol/32 = A(u1,u2;u4,u5) + A(u5,u4;u2,u1) with A = Σ_{m≥0} (−u5)_m/m! · u4/(u4+m) · u2 B(u2+u4+m, 1+u1+u3), from splitting the polytope at t2 = t4 (min-integral), B the beta function. The m=0 term is a Γ-product times the rational constant a2/(a2+a4); all irreducible MZV content is in the m≥1 part.
4. **Reversal** (a1..a5)→(a5..a1) is a symmetry (t2↔t4).

**Checks** (`numcheck_pt1.pkl`, `numcheck_pt2.pkl`, `numcheck_dbl.pkl`):

| check | residual | digits |
|---|---|---|
| t-integral vs direct x-coordinate integral, pt (1/2, 1/3, 1/4, 1/5, 1/6) | 8.1e-52 | 51 |
| same, pt (0.3, 0.8, 0.5, 0.6, 0.4) | 0.0 at dps 80 | ≥72 |
| anchor u=(1,1,1,1,1): Vol/32 = 2/15 (exact hand value ∫(1−t²)²/4) | 0 | exact |
| A+A_rev series vs t-integral, pt 1 | 2.4e-22 | 22 |
| double series (recurrence + double Richardson) vs t-integral, both pts | 8e-13 / 2.3e-12 | ~12 (Richardson-limited) |

### C1. Exact engine (`engine5.py`; EXACT computer algebra, not certified)

Vol/32 = G̃·Σ_{m≥1} α_m + (m=0 Γ-products), G̃ = Γ(1+u1+u3)Γ(1+u2+u4)/Γ(1+u1+u2+u3+u4). The m≥1 Euler sums reduce to MZV words via Σ_m m^{−c} H_{m−1}^{(w)} = ζ(c,w) and stuffle; reduced with CHN1's DS tables (nullities match conjectured values; W=11,12 via flint). Symbolic coefficient forms c_W(a) = [u^W]Vol/32 through **W=10** in `forms10.pkl` (dicts (e1..e5)→Fraction per basis monomial). Two algebra bugs were caught and fixed by the controls (the p_n closed form must subtract a4^{n−1}; q_m = (−1)^m(Y^m − X·Y^{m−1})).

**Controls — ALL PASS:**
- u1=0 ⇒ CHN1's P4 coefficients of (a2,a3,a4,a5) exactly; u5=0 ⇒ P4 of (a1..a4) (5 points, W≤10, full coordinate vectors including products; the m=0 part added at points).
- u3=0 ⇒ [Γ(1+u1)Γ(1+u2)/Γ(1+u1+u2)]·[Γ(1+u4)Γ(1+u5)/Γ(1+u4+u5)]; u2=0 ⇒ Γ(1+u4)Γ(1+u3+u5)/Γ(1+u3+u4+u5); u4=0 ⇒ Γ(1+u2)Γ(1+u1+u3)/Γ(1+u1+u2+u3) (W≤10).
- Reversal invariance: exact polynomial identity of every coordinate form, W≤10.
- **Cauchy/FFT numeric check** (N=64, r=1/30 = R/6 with R=0.9/max|a_i|, dps 70, 64 complex evaluations of the A+A_rev series): c0 = 1; c2..c8 at generic direction (1,2,4,3,5) match the exact engine values with relative differences **1.4e-34, 7.3e-23, 1.0e-21, 6.2e-21, 1.9e-20, 3.3e-20, 1.6e-20** (`fft_check.pkl`).

### C2. The irreducible forms (EXACT)

Basis conventions as CHN1. **Factored results:**

- **Weight 8, ζ(3,5):** P8^(5) = **a2a3a4 (a1+a2+a3+a4+a5) C4 / 10**, C4 an **irreducible quartic** (55 terms; sp.factor). C4(0,a2,a3,a4,a5) = a5·C3(a2,a3,a4,a5) and C4(a1,a2,a3,a4,0) = a1·C3(a1,a2,a3,a4) exactly (P4's smooth cubic).
- **Weight 10, ζ(2)ζ(3,5):** −**a2a3a4 (Σa) (a2+a4)(a1+a3+a5) C4 / 10** — the same C4; generalizes P4's (a2+a4)(a1+a3) to (a2+a4)(a1+a3+a5).
- **Weight 10, ζ(3,7):** **a2a3a4 (Σa) Q7′ / 224**, Q7′ an **irreducible sextic** (182 terms, integral), with Q7′(0,a2..a5) = a5·Q5(a2..a5) and Q7′(a1..a4,0) = a1·Q5(a1..a4) exactly (P4's quintic).

**Divisibility:** a2, a3, a4 and Σa divide all three forms; a1 and a5 do **not** (consistent with the degenerations: u1=0/u5=0 → P4 which has irreducibles; u2/u3/u4=0 → Γ-products). So on positive directions the irreducible-free locus through weight 10 is exactly **{C4 = 0} ∩ {Q7′ = 0}**. Unlike P4 there is no cubic — the P5 irreducible carrier is a quartic (W=8) and a sextic (W=10).

### C3. The irreducible-free locus — PARTIAL (ran out of turns)

- **Negatives (EXACT):** none of the candidate pair-symmetric families lies in the locus: C4 ≢ 0 and Q7′ ≢ 0 on (x,x,y,y,z), (x,y,x,y,z), (x,y,z,y,x), (x,x,y,x,x), (z,y,y,x,x), (z,y,x,y,x), (x,y,y,x,z), (x,x,x,y,z), (x,y,x,x,z). **The diagonal is not irreducible-free for P5:** C4(x,x,x,x,x) = −232x⁴ (in stark contrast to P4, where the diagonal kills C3).
- **Boundary:** on the face a1=0 the locus is P4's Π1∪Π2 in (a2..a5); on a5=0 similarly; on a2/a3/a4 = 0 the volume is a Γ-product for all a (all-orders irreducible-free, from the degeneration formulas).
- **Interior (NUMERICAL):** scanning 6809 positive points of {C4=0} (a3=1 slice, random (a1,a2,a4), quartic roots for a5>0): Q7′ takes **both signs** (1628 positive, 5181 negative), so {C4 = Q7′ = 0} has **nonempty positive (curved) components**; the closest sampled approach was |Q7′|/|a|⁶ ≈ 3e-5. The precise components (resultant/Gröbner analysis) and whether they admit finite Γ-sum closed forms were **not determined** — this is the main open item, together with C4 (the six-point disk-integral comparison).

**Γ-sum loci:** beyond the boundary faces (where P4's (B+G0)/2 theorem applies on Π1∪Π2), no Γ-sum locus was found in the interior; none of the tested linear families admit one. Whether the curved {C4=Q7′=0} components admit a finite Γ-sum is untested.

### Labels
- **EXACT** (computer algebra; engine validated by the controls above but not certified): all coefficient forms, factorizations, restriction identities, divisibility, candidate-family vanishing tests, degeneration controls.
- **PROVED** (elementary derivations given in the text): the volume formulas (1)–(3), the degeneration formulas, reversal symmetry.
- **CONJECTURAL:** MZV independence (needed only to read "coefficient ≠ 0" as "irreducible genuinely present" for real numbers); all DS reductions used only proved relations.
- **NUMERICAL:** the volume-representation checks (digit counts above), the FFT check, the Q7′ sign scan.

### Files and reruns
- `engine5.py` — the P5 symbolic engine; rebuild W≤10 forms: `python -c "import engine5,pickle; f,e=engine5.reduced_forms(10); pickle.dump(f,open('forms10.pkl','wb'))"` (≈5 min for build + reductions; needs `mzv_exact.py`, `kimilib.py`).
- `forms10.pkl` — all DS-basis coordinate forms, W=3..10; `forms8_W8.pkl` (W≤8 checkpoint).
- `kimilib.py` — m=0 Γ-product contribution at points, Γ-product series, evaluation helpers.
- `partB_planecoords.pkl`, `numcheck_pt1.pkl`, `numcheck_pt2.pkl`, `numcheck_dbl.pkl`, `fft_check.pkl` — checkpoints.
- Harness step files: Part B machinery/comparison `step_13.py`–`step_18.py`; volume-formula numerics `step_28.py`–`step_41.py` (tanh-sinh quadrature, A-series, double series); engine build/controls `step_44.py`–`step_64.py`; FFT check `step_61.py`–`step_67.py`; forms analysis `step_69.py`–`step_75.py`.
- `progress.log` has one line per finished step; `summary.md` "Kimi job" section updated.

**Bottom line.** For P5 the weight-8 irreducible enters via P8^(5) = a2a3a4(Σa)·C4/10 with C4 an irreducible quartic restricting to P4's C3 on both boundary faces (times the remaining variable); the weight-10 forms are −a2a3a4(Σa)(a2+a4)(a1+a3+a5)C4/10 (ζ(2)ζ(3,5)) and a2a3a4(Σa)Q7′/224 with Q7′ an irreducible sextic restricting to P4's Q5. The irreducible-free locus in the interior of the positive orthant is {C4=0}∩{Q7′=0}: nonempty and curved, containing none of the natural pair-symmetric linear families — a genuinely different picture from P4's two planes.