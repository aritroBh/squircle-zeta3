# CX9 independent referee report

**Overall verdict: PROOF COMPLETE, with the corrections in item 6.**

The symbolic cache bridge and the independent k=5 boundary certificate are supplied in this review. All requested original symbolic certificates have rerun successfully. The theorem concerns admissible depths `0<=delta<=min(6,2d-4)`, off-diagonal keys from k=d, and stable diagonal keys from k=2p+1. The actual diagonal initial value at k=2p is retained separately; E3 establishes blindness of its adjacent child to that exception. The resulting weight-five cubic holds for every integer k>=5, with constants defined by WT5's exact kernels/integrals. This does not prove WT5's separate numerical classical closed-form relations.

Only this review folder was modified. All reported fresh checks were actually run; their scripts and logs are retained here.

## 1. Fitted e(5,1), e(6,1), e(6,2) and completeness of the Newton system — OK

The discovery of these expressions by fitting is immaterial after verification of the full identity set. The complete new Newton-index set is:

| Depth delta | Every required s | Non-adjacent p>=q+2 | Adjacent (p,p-1), including fold |
|---|---|---|---|
| 5 | 0,1,2 | all symbolic zeros in d,m | all symbolic zeros in p |
| 6 | 0,1,2,3 | all symbolic zeros in d,m | all symbolic zeros in p |

All seven original `prove_V56.py DELTA S` runs pass. All seven adjacent identities also pass using the submitted `eR` and `feed_terms` expressions evaluated in the exact rational-function field `QQ_I.frac_field(p)` (`check_adj_fast.py`, `adj_fast.log`). This changes the simplification method, not the candidate or feed. Original adjacent expansion-script rerun status is recorded in item 3.

There is no omitted Newton index: above s=floor(delta/2), both sides vanish because every parent depth delta'<=delta has e(delta',s)=0. The value, intermediate slope, and top identities are all present. Lower depths use the corrected CX8 theorem.

For diagonals, define the stable coefficients directly by the sole off-diagonal parent:

`ehat_pp(delta,s) = (p+1) sum_{delta'<=delta} tau1_(delta-delta') e_(p+1,p)(delta',s)/(4p-delta'-1)`.

Here e and ehat denote absolute coefficients; the library e is normalized by rho(0). This construction satisfies the diagonal T1 equation identically in the basis `binom(k-2p-1,s)`. It requires no additional equation falsely identifying ehat(delta,0) with the diagonal initial rho(delta). Feeding an adjacent child converts the basis by

`Ediag(delta,s) = sum_{t>=0} (-1)^t ehat(delta,s+t)`.

The submitted adjacent certificates implement both this alternating conversion and the full fold. Their fold powers `(-i)^(-4-t)` and `i^(-6-delta)` follow from removing integer multiples of 4 in the original p-dependent exponents.

Normalization and conjugation matter here. The adjacent script divides by the real rho_(p,p-1)(0), with parent ratios

- `rA0 = i(4p-3)/(p+1)`;
- `rB0 = (4p-1)(4p-3)/(p(p+1))` (two generations).

The E3 script instead divides by the purely imaginary rho_pp(0). Its parent ratio is `2i(4p-1)/(p+1)`, and physical conjugation of the discrepancy gives `Dc=-Delta(conj=True)`. The submitted signs are correct. Fresh checks of tau against its DL definition and of the actual boundary feeds support this audit (item 4).

**Why this proves the fitted forms.** Fix an arbitrary finite k and start at degree k with the true initial data. Descend finitely in degree. Off-diagonal initial values agree by e(delta,0)=rho(delta). At stable off-diagonal keys the Newton identities make the candidate satisfy T1; at stable diagonals the single-parent construction does. At a diagonal initial key keep its actual rho; E3 allows its adjacent child to receive the stable extrapolation instead. The finite triangular recursion has a unique solution, so the candidate equals the actual tower. This is precisely CX8's descending-induction argument, without an infinite-degree recursion or an assumption about uniqueness of a fit.

## 2. Initial rho(5), rho(6), log derivatives and Bell assembly — OK

I read `result.md` and `proof_depth6.md` fully and checked the actual library/certificates. Independently of their L cache, `check_initial.py` uses

`A(v) = C v^(2d-1) S_m(t(v)) (v+i)^(2-2d)/(v^4-6v^2+1)`,

where C is independent of v and `t(v)=v/4+1/2-1/(4v)`. It differentiates the rational logarithmic derivative of the algebraic factor and expands log S(t(i+x)) using t'(i)=0. Put Rj=S_m^(j)(t*)/S_m(t*), ell_r=(log A)^(r)(i). The resulting derivatives are:

```
ell1 = -i(d-2)
ell2 = (3d-3)/2 - i R1/2
ell3 = i(7d-5)/2 + 3 R1/2
ell4 = (-45d+21)/4 + 6i R1 - 3(R2-R1^2)/4
ell5 = -(93i/2)d + 21i/2 - 30 R1 + (15i/2)(R1^2-R2)
ell6 = (945d-45)/4 - 180i R1 + (135/2)(R2-R1^2)
       + (15i/8)(R3-3 R1 R2+2 R1^3)
```

With Lr=ell_r/r!, all six agree **symbolically in d,m** with `L_parsed.pkl`. This repairs the submission's insufficient random-point-only comparison of that cache to the proved construction. The Bell recurrence

`a0=1; n a_n = sum_{r=1}^n r Lr a_(n-r)`

then proves a1..a6 equal the submitted arho formulas. Only R1..R3 are needed through order six; R4 first enters at order eight. The two displayed a5,a6 numerators in proof_depth6.md were separately transcribed and checked against the library, also with zero symbolic residuals (`compare_display.py`, `display.log`).

`check_S_symbolic.py` independently derives R1..R4 by expanding S(t*+h)/S(t*), using `q2(t*+h)=2ih+2h^2` and the Stirling-tail polynomials through order four. It checks their recurrence and anchors, cancels falling-factorial ratios before evaluating endpoints, and obtains all four library Rj as exact rational identities in d,m. Evidence: `S_symbolic.log`, four zero residuals.

Additional fresh exact checks (`initial.log`):

- **116/116** finite-polynomial R1..R4 checks, every canonical key of degrees 1..9, including q=0, adjacent and diagonal shapes.
- **24/24 each** for rho(5),rho(6), against the specified `TOWN1_tower_large_n/data/towers_D8_K52.pkl`: d=5,6,7,8,12,20,40,52, p in {d,d-1,ceil(d/2)}.

These compare the actual initial coefficient with rho(0)*a_delta. The off-diagonal factor two is already in rho(0), and must not be applied again. These are fresh spot checks, not the submission's historical 720-cell totals.

## 3. Genuinely symbolic reruns and domains — OK

The symbolic variables remain unevaluated: d,m for non-adjacent identities, p for adjacent identities and E3. All sums in these finite-depth certificates have fixed finite ranges; none samples a range of degrees as a substitute for proof.

| Certificate | Cases | Fresh result |
|---|---|---|
| Original prove_V56.py | (5,0..2),(6,0..3) | all seven zero numerators; symbolic.log |
| Original prove_adj56.py | (5,0..2),(6,0..3) | all seven zero remainders; adj_5_0.log, adjacent.log |
| Same submitted adjacent feeds in QQ(i)(p) | all seven | all zero, individual runs 0.13..26.42s; adj_fast.log |
| Original prove_E3b56.py | 5,6 | both zero numerators; symbolic.log |

Original adjacent internal simplification times were (5,0): 169.7s; (5,1): 20.7s; (5,2): 1.5s; (6,0): 496.3s; (6,1): 93.9s; (6,2): 11.3s; (6,3): 0.6s. Every original process returned status 0. The two value cases are slow because the original checker expands a product of all term denominators; the independent exact-field evaluations avoid that expansion.

**Domain correction.** The defining bound is delta<=2d-4. Thus **both depths 5 and 6 begin at d=5**, and adjacent d=2p-1 begins at **p=3**. The briefing's d>=6 / p>=4 thresholds are sufficient interior ranges, not the sharp lower bounds. Restricting the conclusion to those ranges would omit the degree-five case essential to W5. The symbolic rational identities have no singularity at d=5 or p=3: the relevant Pi factors and feed denominators are nonzero, as are p and p+1. The new formulas are therefore valid there as well.

Independently, `check_boundary.py` compares all three degree-five keys `(5,0),(4,1),(3,2)` at depths 5,6 and k=5,6,7,8,9,12,21,40,52 with the saved exact towers: **54/54** matches. This is corroboration of the boundary, not an extrapolation premise.

## 4. Blindness E3 and first diagonal levels — OK

Both original E3 certificates pass as identities in p, with zero numerator at delta=5 and delta=6 (internal simplification times 7.3s and 15.4s). They include the direct DL2 contribution and the conjugated fold.

`check_boundary.py` separately implements the T1 diagonal feed using exact Q(sqrt2,i) arithmetic, with tau computed from the DL rational functions rather than the table. The table agrees through n=7. Fresh comparisons against the specified TOWN1 pickle pass:

- **336/336** cells: adjacent and diagonal keys for p=3,4,5,6, depths 0..6, first six levels. At diagonal initial levels the comparison uses the independently certified rho formula; later levels use ehat.
- **8/8** zero direct-plus-fold feeds of the initial-minus-stable-extrapolated discrepancy, p=3..6 and delta=5,6. All eight corresponding diagonal discrepancies themselves are nonzero.

The first diagonal levels are thus treated correctly; E3 asserts blindness, not equality of the exceptional initial point and the stable polynomial. Evidence: `boundary.log`.

The new top stable-diagonal coefficients, normalized by rho_pp(0), simplify symbolically to

```
ehat(5,2)/rho_pp(0) = -4(1+i)(p-2)/((4p-5)(4p-3))
ehat(6,3)/rho_pp(0) = 5(-1+i)/((4p-7)(4p-5)(4p-3))
```

They are nonzero throughout the admissible p>=3 range. The off-diagonal tops `-2i(d-4)/Pi2(d)` and `-5/Pi3(d)` are nonzero for d>=5. Hence the new tower degrees are exactly 2 and 3, not merely upper bounds.

## 5. The k=5 boundary and WT5 identification — OK

This was my own fresh run, independent of the pre-screen and TOW3R. `check_kernels.py` uses only Python Fraction arithmetic on the full saved partial-fraction dictionaries from `kern_5.pkl` through `kern_9.pkl`.

For **each** key `(4,0),(3,1),(2,2)` it verifies the full rational-kernel equality

`K5 = 4 K6 - 6 K7 + 4 K8 - K9`.

There are **40/40 real-residue QSI scalar components**, all zero residuals, including absent-pole coordinates through the union of supports. The +/-i coefficients and constant parts are also compared and are all zero separately at each k=5,6,7,8,9. This checks the entire saved kernel, not only selected real residues.

For every key, the first three differences at k=5 equal WT5 `deltas_deg4.json` exactly:

```
Delta K5   = K6-K5                    = D1  -> integral E1
Delta^2 K5 = K7-2K6+K5                = D2  -> integral E2
Delta^3 K5 = K8-3K7+3K6-K5           = D3  -> integral E3
```

All **nine** difference-kernel comparisons pass (`kernels.log`). WT5's `wt5_B_gexpr.py` identifies those D1,D2,D3 integrals with E1,E2,E3; I inspected that construction without executing its external writes.

**Consequence.** Degree-five towers exhaust depths 0..6. WTS Lemma 3 therefore gives degree-four kernel residues polynomial in k of degree at most three for k>=6. The same bound applies to each +/-i simple-pole coefficient; its four zeros at k=6,7,8,9 force identical vanishing. The fixed real-linear tower-to-kernel map and fixed integration map preserve polynomiality. The exact k=5 identity above puts that level on the same cubic. Its differences have just been identified with WT5, so

`W5(k) = W5'' + (k-5) E1 + binom(k-5,2) E2 + binom(k-5,3) E3`

for **every integer k>=5**, with W5''=W5(5). The exceptional k=4 is outside this conclusion. This audit does not prove any separate integer-relation formula for the constants.

## 6. Wording and reproducibility fixes — OK, corrections required

1. Replace the random-point L-cache comparison in proof_depth6.md section 1 by the symbolic bridge supplied here. Samples alone do not certify equality of that cache with the proved construction.
2. State sharp admissibility and level ranges: delta<=min(6,2d-4), d>=5 for each new depth, adjacent p>=3; off-diagonal k>=d, stable diagonal k>=2p+1, actual initial diagonal k=2p retained separately.
3. Distinguish normalized library a_delta/e from absolute rho/e in T1, and display the ehat shift, parent ratios, and conjugation signs in item 1. Never claim E3 puts the initial diagonal value on the polynomial.
4. Cite every new value, slope, and top equation plus E3 and finite descending induction to certify the fitted forms. Unused depth-7/8 entries in e_fits.pkl remain candidates. Non-adjacent all-s top equations do not prove the all-depth law.
5. Include the explicit k=5 extrapolation, the four imaginary-pole zeros, and the WT5 difference-kernel identification in the all-k consequence. Lemma 3 as stated starts at k=6 for degree-four kernels.
6. Call the tower-to-kernel operator fixed real-linear, or include conjugate coordinates explicitly; it is not complex-linear in the +i coordinates alone. Polynomial degree is preserved either way.
7. Upgrade the cubic structure to a theorem while retaining the independent numerical status of WT5's classical closed-form relations. This review does not upgrade those relations.
8. Paper `prop:w5k` currently states 5<=k<=52. Its structural range can become every integer k>=5, citing TOW3 and this boundary certificate, retaining the exceptional k=4 and the status labels for the constants. No paper file was edited here.
9. The advertised `verify_tow3.py` is absent. Give the actual certificate commands with required DELTA/S arguments. Historical `step_33.py` incorrectly applies the off-diagonal factor twice in its expected ratio; remove its extra `mult` if maintaining it. It was inspected, not executed here, and is not load-bearing. Avoid code-generating historical steps as read-only verification commands.
10. Remove the duplicated Task D heading and reconcile its status text. Finite validations and the conjectural sub-top family do not prove an all-d M_d formula. Task D is outside this verdict.

## Reproducibility and limits

Use `<project>/.venv/bin/python -B`. `run_proofs.py` invokes original submitted certificates by absolute path, with cwd here and bytecode disabled. Without arguments it runs all seven V cases and both E3 cases; with DELTA S it runs the corresponding adjacent case. `run_adj_rest.py` runs the remaining six original adjacent cases; (5,0) was run separately. Logs retain every result and timing. `check_adj_fast.py` provides a much shorter exact-field evaluation of all seven unmodified adjacent feeds. The other local check scripts reproduce the independent S, log/Bell, display, boundary, and kernel checks described above.

I did not regenerate the large TOWN1 pickle, rerun the historical 720/414/38939-cell totals, rerun every corrected CX8 depth<=4 certificate, redo WT5's integration/closed-form relation certificates, or certify Task D. The corrected CX8 theorem and the written WTS tower-to-kernel argument are inherited premises. No dependency was installed, and no source or data outside this folder was modified.
