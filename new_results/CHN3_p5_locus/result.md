# CHN3: the 5-chain P5 — irreducible-MZV vanishing locus

**Work directory:** `<project>/new_results/CHN3_p5_locus/`

## Executive summary

Step 1 (the exact geometry of the weight-≤10 irreducible-free locus **V10 = {C4=0} ∩ {Q7′=0}**) is **complete and rigorous**. Steps 2–4 are **partially done**: a validated exact scalar point-engine was built, the W=11/12 double-shuffle tables were constructed (nullities 9 and 12, the conjectured values), and W=11/12 coefficient values were computed at a point, but the full symbolic W=11/12 forms, the intersection, and the Γ-sum test were **not finished** (the symbolic `build_engine(12)` proved too slow under the machine load of ≈174, and my vectorized modular engine has one identified bug, now understood).

---

## Step 1 — Exact geometry of V10 (COMPLETE)

**Setup.** From CHN2's `forms10.pkl` I extracted (exact polynomial division, verified `a2a3a4(Σa)·C4/10 == P8`):
- **C4**: homogeneous quartic, **55 terms**, irreducible over ℚ, reversal-invariant, C4(1,1,1,1,1) = **−232**; deg 4 in a1,a5 and deg 3 in a2,a3,a4; leading coefficient in a5 is the constant **4**.
- **Q7′**: homogeneous sextic, **182 terms**, irreducible over ℚ, reversal-invariant; deg 6 in a1,a5 and deg 5 in a2,a3,a4; leading coefficient in a5 is the constant **64**.
- Confirmed exactly: `P10(ζ(2)ζ(3,5)) = −a2a3a4(Σa)(a2+a4)(a1+a3+a5)·C4/10` (so the weight-10 ζ(2)ζ(3,5)-form vanishes wherever C4=0 — no extra condition).

Since both leading coefficients in a5 are **constants**, the a5-resultant loses nothing. Dehomogenising a3=1 (every positive point has a3>0, so nothing positive is lost):

**Result (EXACT):**
$$\operatorname{Res}_{a5}(C4,Q7') = 16\cdot a1\cdot(a1+a3)\cdot G(a1,a2,a3,a4),$$
with **G irreducible over ℚ of degree 22** (2138 terms dehomogenised). Files: `R_a1a2a4.txt`, `R_factors.txt`.

**Projection of V10** to (a1,a2,a3,a4)-space is `{a1=0} ∪ {a1+a3=0} ∪ {G=0}`:
- `{a1=0}`: boundary face (carries the stratum {a1=0,a5=0} and the P4 curve {a1=0,C3=Q5=0}) — **not** in the open positive orthant.
- `{a1+a3=0}`: **no positive points** (a1+a3>0 on the orthant).
- `{G=0}`: the **only** component meeting the open positive orthant.

**Birational model (EXACT):** the degree-1 subresultant gives, on {G=0} with A≠0,
$$a5 = -B(a1,a2,a4)/A(a1,a2,a4),$$
A (692 terms), B (748 terms) in (a1,a2,a4); the degree-0 subresultant equals R exactly (`subres0.txt == R_a1a2a4.txt`, verified). So **V10's positive part is birational to the single irreducible degree-22 surface {G=0}** — one ℚ-irreducible component (projectively 2-dimensional).

**Positive point (EXACT certificate):** the point
$$(a1,a2,a3,a4,a5) = (1,\,1,\,1,\,a4^\ast,\,a5^\ast),\quad a5^\ast=-B(1,1,a4^\ast)/A(1,1,a4^\ast),$$
where a4\* = unique real root of G(1,1,a4) in (1,2) (isolated to 2⁻⁶⁶ by `dup_refine_real_root`), is a **certified positive point of V10**: exact interval arithmetic (Fraction) gave
- a4\* ∈ (1.293668784620243, …) > 0,
- A(1,1,a4\*) ≈ −2.03·10¹⁵ ≠ 0 (interval excludes 0),
- a5\* ∈ (0.07051999939180706, …) > 0.

Second positive point: a4\* ∈ (5,6) (≈5.1967), a5\* ≈ 6.5304 (numeric). Non-positive lifts: a4=1 (double root) → a5∈{0,−1}; a4≈3.4362 → a5≈−2.083.

**Conclusion of Step 1:** V10 ∩ {open positive orthant} is **nonempty**, lies on a **single ℚ-irreducible component** (birational to {G=0}), and is **projectively 2-dimensional** (a surface, not a curve). This is the locus through weight 10; whether it survives W=11/12 is the subject of Step 3 (not completed).

Labels: **EXACT** (sympy resultant/factor/subresultants; the engine forms themselves not certified) + **NUMERICAL** (mpmath root values).

---

## Step 2 — engine to W=11/12 (PARTIAL)

**Done:**
- `point_engine.py`: exact scalar point-engine (Fraction) mirroring engine5. **Validated EXACTLY against `forms10.pkl` at W≤10 (6 rational points, max |diff| = 0).** One self-corrected bug (a duplicated accumulation line that doubled the regular part).
- W=11 and W=12 DS tables build correctly: **nullity 9 (W=11)** and **nullity 12 (W=12)**, matching the conjectured MZV dimensions; `DSTableFlint` used for W≥11.
- Computed c11, c12 at the test point (1,2,3,4,5), saved in `test_pt_w1112.pkl`:
  - c11 = [276894734985, −84576432720, −23142001248, −69057414816/7, −181411559208/35, −1863957600, 17474375160, −218116800, 109753920] (basis [ζ11, ζ2ζ9, ζ2²ζ7, ζ2³ζ5, ζ2⁴ζ3, ζ2ζ3³, ζ3²ζ5, **ζ3ζ(3,5)**, **ζ(3,3,5)**]); so at this point ζ(3)ζ(3,5) has coeff **−218116800** and ζ(3,3,5) has coeff **109753920**.
  - c12 irreducible coords: ζ(2)²ζ(3,5)=552318912, ζ(2)ζ(3,7)=171509400, ζ(3,9)=1022902694352/5, ζ(5,7)=−654673492512/5, ζ(2,1,1,8)=1722729600.
- `modengine.py`: vectorized modular engine (numpy int64 mod p, K points per pass). Validates at W=3; **bug identified**: my `st_dicts` drops the stuffle multiplicity `c` (should be `c1*c2*c`); this first bites at W≥5 (e.g. the word (3,1,1) came out −15 instead of −9). Fix is a one-line change (`nd[w]=(nd[w]+c1*c2*c)%p`), **not yet applied/validated**.

**Not done:** the full symbolic W=11/12 forms and their factorisations. The symbolic `engine5.build_engine(12)` (pure-Python Fraction stuffle arithmetic) did not finish within ~26 min under system load ≈174 and was abandoned in favour of the point/modular engines. The controls (faces a1=0/a5=0 vs CHN1's P4 forms, a2/a3/a4=0 → Γ-products, reversal, FFT of c11/c12) are scripted (`fft_check_1112.py`) but not run.

## Steps 3–4 — NOT done
The intersection of V10 with the W=11/12 obstruction forms, and the Γ-sum/Stokes analysis, remain open. The planned method (exact, no Gröbner): test ideal membership of each obstruction form F in the complete intersection (C4,Q7′) by solving F = U·C4 + V·Q7′ as Fraction linear algebra (`linalg_tools.py` written, untested); on {G=0}, F restricts via a5=−B/A.

---

## P5 theorem candidate (status)

Through weight 10 (**EXACT**, modulo engine certification): for positive a, the volume of the 5-chain contains an irreducible MZV (ζ(3,5) at weight 8, or ζ(3,7) at weight 10) **unless** (a1,…,a5) lies on V10 = {C4=0}∩{Q7′=0}, which is a nonempty irreducible surface in the positive orthant (positive points proved to exist). Whether irreducibles reappear at weight 11/12 on V10 — i.e. whether V10 survives to all orders as an irreducible-free locus — is **undetermined**. The full "P5 theorem" (analogue of P4's "Vol is an irreducible-free Γ-sum iff a ∈ Π1∪Π2") is **not established** here.

## Labels
- **EXACT** (computer algebra, not certified): all of Step 1 (resultant, factorisation, subresultants, certified positive point via exact interval arithmetic); the scalar point-engine validation (max diff 0); the DS-table nullities 9 and 12.
- **PROVED** (from CHN2, reused): the volume representations (t-integral, h(t), A+A_rev).
- **CONJECTURAL**: MZV independence (needed only to read "coefficient ≠ 0" as "irreducible genuinely present" for real numbers); all DS reductions used only proved relations.
- **NUMERICAL**: the mpmath root values of a4\*, a5\*; the W=11/12 point values pending FFT cross-check.

## Evidence table
| check | method | result |
|---|---|---|
| P8 = a2a3a4(Σa)C4/10 | sympy exact division | exact |
| P10(ζ2ζ(3,5)) = −a2a3a4(Σa)(a2+a4)(a1+a3+a5)C4/10 | sympy | exact |
| C4, Q7′ irreducible / reversal-invariant | sympy factor / simultaneous subs | exact |
| Res_a5 = 16·a1·(a1+a3)·G, G irred deg 22 | sympy resultant (49.5s) + factor (5.2s) | exact |
| subres0 == R | sympy | exact |
| a5 = −B/A reproduces both positive points | nroots vs subresultant | residuals 2.4e-35 / 4.6e-32 |
| positive point (1,1,1,a4\*,a5\*) | dup_refine + Fraction interval arithmetic | a4\*>0, a5\*>0 **certified** |
| scalar point-engine vs forms10, W≤10, 6 pts | exact Fraction comparison | **max diff 0** |
| DS tables W=11,12 | flint rref | nullities 9, 12 |

## Files and reruns
- `common.py` (C4/Q7′ extraction), `C4.txt`, `Q7p.txt`.
- `run_resultant.py` → `R_a1a2a4.txt`, `R_factors.txt` (rerun: `nice -n 10 PY run_resultant.py`, ≈55 s).
- `run_subres.py` → `subres1_A.txt`, `subres1_B.txt`, `subres0.txt` (≈47 s).
- `point_engine.py` (validated scalar engine), `modengine.py` (modular engine; fix `st_dicts` to include `*c`).
- `linalg_tools.py` (interpolation/ideal-membership, untested), `fft_check_1112.py` (unrun).
- `test_pt_w1112.pkl` (c11/c12 at (1,2,3,4,5)), `progress.log`, `summary.md`.
- To resume Step 2: fix `st_dicts` in `modengine.py`, re-validate vs forms10 at W≤10, then compute W=11/12 coords at enough points to interpolate the quotient forms (after dividing by a2a3a4(Σa)) and run the factor/membership tests.