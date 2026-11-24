# CX4 independent adversarial review

**Verdict: GAP in the submitted proof, at Lemma B's adjacent branch.** Its displayed recursion is false: it counts the diagonal feed twice. `step_39.py` repairs this, but `result.md` retains the error and misreports `step_38.py` as passing. The correction below supports an all-k proof; this is not a counterexample to Conjecture 3. I did not consult the first reviewer's files.

## Lemma A — valid on the stated reducer's domain

Fix k and induct downward in degree d. The initial degree-k coefficients have order at most 2k−2, by the square-root expansion below. For a parent of degree d+1, its higher part has order at most 2d. Crucially, **Int receives only rational principal parts of order r≥2**, so
\(\int c(v-i)^{-r}dv=-c(v-i)^{1-r}/(r-1)\).
There are no logarithmic terms: simple poles and the log monomials are excluded from Int. The unrestricted assertion about integration would be false.

Directly, \(DL_1=(v^2+1)/[v(v-b)(v+a)]\), with \(DL_2\) obtained using \((v-a)(v+b)\). Both have simple zeros at ±i; their derivatives at +i are \(\tau_1=(-1-i)/2\), \(\tau_2=(-1+i)/2\). Evaluation at b adds only a constant. The map 1/v exchanges ±i with nonzero derivative, so F preserves their paired pole-order bound. Thus each feed lowers the bound twice, proving A for d≥2. At d=1, HigherPart and Gm vanish. `step_38` incorrectly treats an empty pole list as order zero when comparing it with −1.

## Lemma B — denominator correct; adjacent numerator wrong

A parent top term has order 2d; integration contributes −1/(2d−1), multiplication by DL contributes its first derivative, and the IBP minus sign cancels. Lower terms cannot reach order 2d−2. At +i, F's leading multiplier is −1 on the coefficient formerly at −i: \((1/v)'|_i=1\), \(i^{-2}=-1\).

Consequently the correct feed numerator is
\[
(i+1)\tau_1h_{i+1,j}
+[i\ge j+2](j+1)\tau_2h_{i,j+1}
+[j=i-1]i\tau_2(h_{ii}-\overline{h_{ii}}).
\]
Divide by 2d−1 and add \(\delta_{k,d}\rho_{ij}\). Examples: (2,0) gives \((3\tau_1h_{30}+\tau_2h_{21})/3\); (1,1) gives \(2\tau_1h_{21}/3\); (2,1) gives \((3\tau_1h_{31}+2\tau_2(h_{22}-\bar h_{22}))/5\).

**Exact counterexample:** at k=4, \(h_{31}=-5\sqrt2i/4\), \(h_{22}=15\sqrt2i/16\), and actual \(h_{21}=-3\sqrt2/4\). The printed B predicts \(-15\sqrt2/16-3\sqrt2i/16\). Actual minus predicted is **\(3\sqrt2(1+i)/16\ne0\)**. My enumeration found 25 failures among 180 cases; the corrected recursion had zero failures.

## Lemma C — all-d derivation works

Put \(T=t(v)\), \(z=(v^2+1)/(2\sqrt2v)\). Expanding the master square root gives
\[
\widetilde R_{d,m}=\frac1{2\sqrt2vT(1-T)}
\sum_{j=1}^d\binom{1/2}{j}z^{2-2j}U_j^{(m)}(T).
\]
This is a formal coefficient identity for every d: each exponential increment starts at degree one. Only j=d reaches order 2d−2. Also \(z'(i)=1/\sqrt2\), hence its indicated power is \(2^{d-1}\).

In \(U_d^{(m)}\), the support conditions \(s\le m\), \(d-s\le d-m\) force **s=m**. Thus \(U_d^{(m)}(T)=2^d\binom dm(1-T)^{2m}T^{2(d-m)}\), giving the claimed value at \(T=(1+i)/2\). The suspicious finite-difference identity is also correct for 0≤m≤d: the polynomial \(s^m(d-s)^{d-m}\) has degree d and leading coefficient \((-1)^{d-m}\); its signed binomial sum is \((-1)^dd!\) times that coefficient.

Finally \(T(1/v)=1-T(v)\), \(z(1/v)=z(v)\), and the Jacobians agree after folding. Swapping the log indices gives identical coefficients, doubling precisely the off-diagonal keys. Using \(\binom{1/2}d=(-1)^{d-1}(2d-3)!!/(2^dd!)\) yields exactly the stated K_d and rho for every d≥2.

## Lemma D — general proof after correcting B

The supplied identity scripts check only finite ranges. General proofs are short:

- **I1:** taking the ratio of the explicit K's gives \(K_{d+1}/K_d=-i(2d-1)/(d+1)\).
- **I2:** both weighted binomials equal \((d+1)\binom dp\); combine I1 with \(\tau_2-\tau_1=i\).
- **I3:** \(\rho_{pp}\) is purely imaginary; \(\rho_{p+1,p-1}=-2p\rho_{pp}/(p+1)\) and \(2p\,i\rho_{pp}=(4p-3)\rho_{p,p-1}\). Substitution proves the adjacent identity for every p≥2.
- **I4:** I1 gives \(\hat h_p=2i\tau_1\rho_{pp}=(1-i)\rho_{pp}\), preserving its imaginary part for every p≥1.

For each arbitrary fixed k, descend from d=k, where the data equal rho. I2 handles non-adjacent keys; the diagonal uses its off-diagonal parent; I3–I4 handle the adjacent key whether its diagonal parent is initial or stabilized. This proves the corrected D for all k, independently of the finite tests.

## Corollary and weight-two identification

At degree two, HigherPart consists only of the double-pole coefficients and their conjugates. The corrected D gives exactly the stated g20/g11 generating series after rational integration and normalization at b. IBP has no endpoint remainder: Gm=O(v−b) at b, and log(t)^i kills every power of log(1−t) at a for canonical i≥1. Simple degree-d kernels have reduction weight d+1, so degree one supplies the weight-two piece in the paper's specified sense.

I independently obtained **zero** for \((1-\epsilon)(W_{10}-\epsilon R_{10})-\epsilon^2Rw2\). This proves the coefficient identity for every k≥2, conditional only on the corrected recursion. Direct exact integration of Rw2·log(t), using the appendix's dilogarithm rules, returned \(\pi^2/24+\log^2(2)/4\), residual **0**. No converse implication is needed.

## Reruns and finite-only evidence

All runs used the project venv Python with `-B`, this cwd, and 900-second alarms; no source/data files were changed. Steps 37, 39, 41, 42 passed their finite checks; 38 exposed the errors above. Step 40 confirmed its coefficient checks, then raised SymPy `CoercionFailed` at its final `apart` call. Step 43's computations, omitting its file-writing log, passed at k=13,14 (71 seconds total).

The actual script ranges are: rho, d=2…12; integrand/fold formulas, d=2…6; U-values, d=2…8; I1–I2, d=2…15 (child degree at most 16); I3, p=2…9; I4, p=1…9; recorded stabilization/kernel checks, k≤12. They do not by themselves prove unrestricted statements.

My printed residuals were **0** for freshly rebuilt initial data/kernels k=2…8, all final keys against the original reducer at k=2,3, finite differences through d=30, and U-values through d=20. These remain finite checks; the unrestricted arguments are given above. A minimal runnable finite-difference check, using the project Python:

```python
from math import comb, factorial
assert all(sum((-1)**s*comb(d,s)*s**m*(d-s)**(d-m) for s in range(d+1)) == factorial(d)*(-1)**m for d in range(31) for m in range(d+1))
```
