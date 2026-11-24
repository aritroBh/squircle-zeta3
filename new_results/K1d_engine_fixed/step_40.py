result = r'''# K1d — Fixed engine predicts the non-local ζ(3) per edge

**Bottom line.** The 1/8→1/4 engine fix is applied and everything the planner predicted is confirmed:
b=0 gives N_e = 159√5/500 per √5 arc (and 0 on the κ=0 arc); the b=1 fresh PSLQ test agrees with the
engine on ALL arcs (cz3 = Z_e + N_e, residuals ~1e-100); and the general single-ratio law is
N_e = Φ_A·ζ₃[A(1−r)] + Φ_L·ζ₃[Li₃(1−r)] with Φ_A, Φ_L explicit closed forms (verified as exact
polynomial identities against the rebuilt pickle, and numerically on all four non-local arcs).

## 0. The bug fix (PROVED)
`K1_nonlocal_law/step_33.py` line 36: the d³ coefficient of √(g0+g1d+g2d²+g3d³) is
`g3/(2Ls) − g1·g2/(4Ls³) + g1³/(16Ls⁵)` (Ls=√g0). Sympy series check: `(8g0²g3−4g0g1g2+g1³)/(16g0^{5/2})`
matches the 1/4 form, not the shipped 1/8 form (`step_40` preamble). Rebuilt → `E3_fixed.pkl`
(40 non-local + 40 local monomials, `step_40.py`).
Check at t=0.3, b=0 edge 1, unweighted logs (`step_41.py`): fixed pickle, first-principles series
engine, and verifier t5 symbolic code all give **E3 = −0.64720423662744655873985211888**
(diffs 3.3e-31, the truncation of the quoted target); old buggy pickle gives 0.5394241874857276 (diff 1.19).

## 1. b=0 (Task 1) — step_42.py
Engine (E3_fixed.pkl, unweighted logs, J-table atoms at x=−κ) vs verifier t3_b0_150.json:

| edge | r | E_arc closed form | cz3 | Z_e | N_e | |diff| vs P3_arc | rel. resid |
|---|---|---|---|---|---|---|---|
| 1 | 2 | √5(13l2³/125 − 72l2²/25 + 333l2z2/250 + 24l2/5 − 24z2/25 + 231z3/500) | 231√5/500 | 18√5/125 | **159√5/500** | 7.45e-102 | 0 |
| 2 | 1/2 | (identical closed form) | 231√5/500 | 18√5/125 | **159√5/500** | 7.45e-102 | 0 |
| 3 | 1 (κ=0) | √2(8l2³ − 24l2² + 3l2z2/16 + 12l2 − 3z2/16 − 15z3/128) | −15√2/128 | −15√2/128 | **0** | 1.25e-102 | 0 |

Per-arc N_e = 159√5/500 on both √5 arcs (r=2 and r=1/2 each); κ=0 arc N_e=0. Verifier PSLQ relations
reproduce with EXACT residual 0 (engine closed form inserted symbolically).

## 2. b=1 (Task 2, fresh prediction) — step_42.py + step_43.py
Engine prediction, then INDEPENDENT PSLQ of X = P3_arc/P0 on [X, ζ3, π²log2, log³2, π², log²2, log2, 1]:

| edge | r | PSLQ relation found (resid) | engine-exact rel. | cz3 = Z_e+N_e | N_e | agree |
|---|---|---|---|---|---|---|
| 1 | 2 | [8788,−615,−297,4748,208,16848,−40560,0] (1.4e-100) | identical | 99√13/4394 + 417√13/8788 = 615√13/8788 | **417√13/8788** | YES |
| 2 | 1/2 | [2000,−267,−132,3336,80,−2160,−7200,0] (3.3e-101) | identical | 99√10/2000 + 21√10/250 = 267√10/2000 | **21√10/250** | YES |
| 3 | 1 (κ=0) | [800,15,−4,−6400,4,19200,−9600,0] (1.0e-102) | identical | −3√5/160 + 0 | **0** | YES |

Agreement on all arcs ⇒ the law predicts (NUMERICAL, ~100 digits).

## 3. The general non-local law for single-ratio edges (Task 3) — step_45..47.py
Edge y(t)=V0+t(V1−V0); term values A0·t, A1(1−t), B(1+κt); κ=r−1; A0=A1−Bκ, A1+B=1;
L2=|V1−V0|², q0=V0·(V1−V0), V02=|V0|²; x = −κ = 1−r. Structural fact (exact on all 6 edges):
the ζ(3)-coefficient of the LOCAL part 6·Lp equals Z_e, so N_e is exactly the ζ(3)-content of 6·Np.
Only the atoms A(x)=ΣHₙxⁿ/n² and Li₃(x) carry ζ(3) in the J-table, giving

**N_e = Φ_A · ζ₃[A(1−r)] + Φ_L · ζ₃[Li₃(1−r)]**  with
  Φ_A = 6B(q0²−L2·V02)/L2^{5/2} · ( 3A0κq0 − L2(2A0 + rA1) )
  Φ_L = 18B(q0²−L2·V02)/L2^{5/2} · ( κ(−A1² + A1B(2κ−1) + 2B²κ)q0 + L2(A0B + rA1²) )
(note q0²−L2·V02 = −L2·h²). ζ₃-contents (re-verified 80 digits, `step_48` inline): A(−1)=−5ζ3/8,
Li₃(−1)=−3ζ3/4 (r=2); A(1/2)=ζ3−π²log2/12 (content 1), Li₃(1/2)=7ζ3/8−π²log2/12+log³2/6 (content 7/8).
At r=1 (κ=0): A(0)=Li₃(0)=0 ⇒ N_e=0. Φ_A, Φ_L were derived by extracting the AA- and L3-coefficients
of 6·Σ_nl coef·J(x) symbolically (`Phi_syms.pkl`), substituting x=−κ, factoring, and simplifying with
the edge constraints; equality with the full pickle is an exact polynomial identity (diff 0, `step_47`).
Checks: Φ_A·(−5/8)+Φ_L·(−3/4) = 159√5/500 (b=0 e1), Φ_A·1+Φ_L·(7/8) = 159√5/500 (b=0 e2),
417√13/8788 (b=1 e1), 21√10/250 (b=1 e2) — all match `step_42` exactly.

## 4. Octagon (two ratios per edge: r=√2 and 1/√2) — what atoms are needed
κ1=√2−1, κ2=1/√2−1; x1=1−√2≈−0.41421356, x2=1−1/√2≈0.29289322; x2=x1/(x1−1) (Landen, diff 3.9e-62),
and 1−x1=√2=1/(1−x2) so −log(1−x2)=+log(1−x1)·(−1) (logs degenerate). The E3 monomials now carry
W1^{c1}W2^{c2}; the new mixed integrals are J(n,a,b,c1,c2)=∫₀¹tⁿ(log t)ᵃ(log(1−t))ᵇlog(1+κ1t)^{c1}log(1+κ2t)^{c2}dt
with (c1,c2)∈{(1,1),(2,1),(1,2)} (pure parts reuse the single-ratio table at x1 and x2).
- weight 2 (1,1): the two-variable dilogarithms **Li_{1,1}(x2,x1/x2) and Li_{1,1}(x1,x2/x1)**.
  Demonstrated exactly (60 digits, inline): with T=Li_{1,1}(x2,x1/x2), T′=Li_{1,1}(x1,x2/x1),
  l_i=−log(1−x_i):  I111 = T(1+1/x1) + T′(1+1/x2) − l1(l2−x2)/x2 − l2(l1−x1)/x1 − (l1−x1)/x1 − (l2−x2)/x2.
- weight 3 (2,1),(1,2): the same Stirling/harmonic partial-fraction machinery produces
  **Li_{2,1} and Li_{1,2} at (x1,x2/x1) and (x2,x1/x2)** plus the two-variable harmonic sums
  H1(x,y)=Σ_{n>m≥1} xⁿyᵐH_m/(nm) (the A-atom analogue), plus one-variable atoms A, Li₃, Li₂ at
  x1, x2, x1·x2 and log products.
- At the octagon point the arguments are Landen-closed, so three-term identities collapse parts:
  observed exactly (PSLQ residual 0.0): Li₂(x1)+Li₂(x2)+l1²/2 = 0 and l2 = −l1 — i.e. the atom basis
  at (x1,x2) is already degenerate; the full weight-3 identification of I21/I12 was NOT completed
  (PSLQ returned only basis-internal identities — the special point collapses the basis).

## Files
step_40.py (rebuild E3_fixed.pkl), step_41.py (t=0.3 check), step_42.py (step_04 rerun → engine_arcs_fixed.json),
step_43.py (fresh b=1 PSLQ), step_45/46/47.py (Φ extraction, factoring, law check), Phi_syms.pkl, Phi_kappa.pkl.
'''
open('result.md','w').write(result)
with open('progress.log','a') as f:
    f.write("step_49: octagon atoms: I111 reduction to Li_{1,1}(x2,x1/x2),Li_{1,1}(x1,x2/x1) verified 60 digits; weight-3 families identified (Li_{2,1},Li_{1,2},H1); octagon-point Landen degeneracies found (Li2(x1)+Li2(x2)+l1^2/2=0). result.md written. DONE.\n")
print("result.md written,", len(result), "chars; progress.log finalized.")
