# Independent adversarial referee report

- **Theorem 12: PROVED AFTER REPAIRS.**
- **Theorem 13: PROVED AFTER REPAIRS.**

The supplied proof contains a false window inequality. A local repair preserves every downstream constant. I found no numerical violation of either theorem. Supplied files were left unchanged; the repairs below are recommendations, not edits to those files.

## Numbered findings

1. **Must fix — `sec_g5_varpi1.tex:179-187`: false sixth-order window bound.**

   At `N=10000`, `k=9999`, `sigma=1/4`, with `q=cosh(sigma)-1-sigma^2/2`:

   ```text
   q-sigma^4/24                 = 3.3946290650949262875085371188614986e-7
   (1.0011/720)*sigma^6         = 3.3945719401041666666666666666666667e-7
   left side minus upper bound = 5.7124990759620841870452194831932445e-12
   required multiplier         = 1.0011168468452749013017176987976823
   ```

   Eight positive points on `sigma=i/4000`, `i=1,...,1000`, violate the stated inequality. The text's geometric argument gives `1/(1-1/896)=1.00111731843575419`, which also exceeds `1.0011`.

   **Repair:** replace `1.0011` by `1.00112`. Then `sqrt(2)*1.00112/720=0.00196638539105965683<0.00197`. Consequently the existing `c_w'`, `T_b`, `E_II`, and both theorem constants remain valid without enlargement. Evidence: `referee_window_violation.out`.

2. **Must fix, numerical evidence — `sec_g5_varpi1.tex:275-277`; `g5_checks.py:17-23,30-35`: inaccurate large-N quadrature.**

   I reproduced the supplied unscaled 40-digit computation, then integrated `exp(theta)*J` using `t=sigma*sqrt(theta)`:

   | N | Quantity | Supplied method, reproduced | Scaled computation |
   |---:|---|---:|---:|
   | 3000 | N^2 exp(theta) abs(J-T) | 1.59577608099707 | 1.59312092855007 |
   | 3000 | N^2 exp(theta) abs(Im(J-T)) | 1.35352106615205 | 1.35024467940650 |
   | 10000 | N^2 exp(theta) abs(J-T) | 1.30957383830269 | 1.57753457852783 |
   | 10000 | N^2 exp(theta) abs(Im(J-T)) | 0.503082404718141 | 0.762246588979991 |
   | 10000 | N exp(theta) abs(J) | 1.78749386806105 | 1.78751705724475 |

   **Repair:** scale before quadrature and replace the table entries. The scaled integrals at 50 and 90 digits differ by `4.33291105071374e-53` at N=3000 and `2.27819226651839e-53` at N=10000. The original quadrature loses accuracy on exponentially tiny integrands. All corrected residuals remain below the analytic bounds; theorem constants are unchanged. Evidence: `referee_unscaled.out`, `referee_grid.out`.

3. **Minor — `sec_g5_varpi1.tex:154`: false assertion on the infinite strip.**

   `abs(epsilon)<1` does not hold throughout `abs(Im sigma)<pi/4`. At N=10000 and sigma=-5, `abs(epsilon)=1.86008310466913879`.

   **Repair:** justify holomorphy using `arg epsilon in (-pi/2,0)`, hence `Re(1+epsilon)>1`. The principal logarithm is holomorphic there. The power series in Step 2 is restricted to the finite rectangle, where its stated smallness bound is valid. No constants change.

4. **Minor — `sec_g5_decay.tex:80`: incorrect global majorant for the logit integrand.**

   For epsilon=-2 and x=1, its absolute value divided by exp(-x) is `1.33384916264797652`, exceeding one.

   **Repair:** for arbitrary real epsilon<1, use `max(1,2^(-1-epsilon))*exp(-x)`. Absolute convergence follows. The subsequent differentiation argument only uses epsilon in (-1/4,1/4), where the original simpler majorant is valid. No constants change.

5. **Minor — `sec_g5_decay.tex:261`: endpoint error in the derivative majorant.**

   The displayed replacement `t^(-1/2)<=0.1` fails for 99<=t<100.

   **Repair:** use `1/sqrt(99)=0.100503781525921208`. Then the derivative coefficient is at most `1/sqrt(99)+sqrt(pi)=1.87295763243143723<1.9`. The later bound by 1.9, and the explicit big-O constant, survive unchanged.

6. **Minor — `g5_consts.py:12-13`: printed endpoint pairs are not rigorous enclosures.**

   Both endpoints are rounded to the nearest displayed decimal. For B_1 the printed pair is `[70.098585335111623752,70.098585335111623752]`, but Arb certifies B_1 is below that displayed lower endpoint by `4.67208797685577282e-19`.

   **Repair:** print the Arb ball with its radius (`v.str(d)`), or round decimal endpoints outward. The internal Arb comparisons are valid; theorem constants are unaffected. Evidence: `referee_interval_display.out`.

7. **Minor — `sec_g5_varpi1.tex:244,293`; `g5_consts.py:89,167`: machine-step/prose mismatches.**

   The script asserts `E_II<71000`, not the advertised `E_II<8000`. Add the stronger assertion; the recomputed value is below 8000. Also, the c_I output label says it uses c_1, whereas the code uses C1=5.1, as intended by the mathematical text. Using the actual c_1 would give `70103.9255399782393`, rather than the implemented `70217.6863816658770`. Correct the label. Neither correction changes the theorem constants.

## Recomputed constants and closure

Independent transcription of the text, using mpmath at 70 digits, is in `referee_checks.py`; printed results are in `referee_checks.out`.

| Quantity | Recomputed value |
|---|---:|
| B_1 | 70.0985853351116237515328 |
| B_2 | 170.657769858025927890213 |
| C_star | 10578907.8235293971071052 |
| (3/(e delta))^3, delta=2-sqrt(pi) | 114.096375458588480449897 |
| c_rest | 114536.103195373587314492 |
| c_1 | 5.088623915831236231239 |
| c_I, using 5.1 | 70217.6863816658769868155 |
| E_II(10000) | 7912.19071213404531373487 |
| Theorem 12 assembled constant | 1209428681.35418576537963 |
| Margin below 1.21e9 | 571318.64581423462037236 |
| Explicit big-O constant (1.21e9+4)/sqrt(pi) | 682669398.349543441398201 |

The six normalized rest pieces at N=100 are:

```text
31918.9426110510885
71215.5143103077935
11125.4793616729538
276.165236285481274
0.00167605595320075035
3.17046096305894402e-10
```

The four E_II terms at N=10000 are:

```text
198.247325187094012
3381.88255343600623
4144.64328477788437
187.417548733060704
```

Strip/model constants independently recomputed:

```text
C_h       = 59.3034578589354477
kappa_K   = 5.89168284930602910
E_1       = 0.0725503508214343118
delta_A   = 1.93122165058864018
e_eta     = 0.00313876533269608724
B_Delta   = 0.0111849075836763186
s_+       = 33.4328122517477757
s_-       = 33.3717598613384232
```

Thus Theorem 13 closes as `71000+120000=191000<1000668`. Theorem 12 closes at the total printed above. The repaired window inequality does not alter either calculation. I reran the supplied Arb and exact SymPy scripts successfully; the exact Lambda_N and varpi_1 identities reduce to zero as claimed.

## Analytic audit and monotonicity

- **Region II, Steps 2-3:** the strip bounds and Cauchy factors are valid. The radius-1/4 circles lie inside the holomorphy strip abs(Im sigma)<pi/4, including a neighborhood of their boundaries. Contact with the boundary of the smaller bounding rectangle is harmless. The real-variable Taylor remainder uses the fourth derivative bound, giving 256*m; no geometric-series remainder factor is missing. Independent 21-by-21 rectangle grids at N=10000, 100000, and 1000000 verified the eta, g, and A bounds. The largest tested Taylor remainder/bound ratio was `0.00662887286946765`.
- **Step 4:** the remainder decomposition, odd-term cancellations, Gaussian moments, and T_a/T_b/T_c factors check after finding 1's repair. The fourth-order coefficient is valid: `sqrt(2)*1.0021/24=0.0590493087855866<0.0591`.
- **Steps 5-7:** model-tail factors, tangent-line tails, the eta expansion and its products, and identification of varpi_1 check. No additional lost factor or sign was found.
- **Modulus lemma:** the D-function derivative proves the stated two-sided bound; log-sum-exp convexity proves convexity of mu. The derivative estimates use secant slopes in the correct directions, and the tangent inequalities justify both tails. A grid with 100 logarithmic N values from 100 to 1e12 and 481 sigma values from -12 to 12 found no failure beyond floating-point roundoff in cancellation-prone differences.
- **Region I:** the exact integral is `2*K_(3/2)(theta)=sqrt(2*pi/theta)*(1+1/theta)*exp(-theta)`. Its prefactor and the negative-tail integration check. Independently scaled quadrature gives `0.628982630402619877` at N=100 and `0.189341504687370159` at N=10000, agreeing with the formula.
- **M_k representation and kernels:** the logit substitution, indentation orientation and residue sign, vanishing vertical segment, differentiation under both integrals, and main-term contour deformation check. Numerical kernel-identity residuals at z=0, 0.3, and 2 were at most `5.44513836812850e-61`.
- **Correction-kernel rotation:** no poles or logarithmic branch cuts are crossed. In the rotated sector, Re L>0 away from the origin; 1+exp(-L) has positive real part there. At L=0 it equals 2. The arc estimates decay exponentially, the ray Jacobian has modulus one, and the factorial cancellation is correct. Checks at k=1,5,99 and r=0,2,20 satisfy the W bound.
- **Integral lemma F:** the geometric lower bound on abs(s+q), logarithmic bound, AM-GM split, and factor 4/lambda are valid, including j=0. A 40-point grid over j=0,1,5,99,400; q=1,2,3,100; lambda=a,3a gave maximum actual/upper-bound ratio `0.0739635759948695`.
- **Rest bounds:** the finite-sum split, factorial tail, D_k and U_3 estimates, and both rest forms check. The four refined-rest logarithmic derivative brackets are negative at N=100 and decrease thereafter for the actual four (c,m) pairs. The exponential-tail pieces also decrease.
- **Every monotonicity reduction:** region I's positive term decreases immediately; dividing its negative-tail denominator by theta proves the asserted ratio decreases. In region II, B_Delta/y^2 is nondecreasing in y; this follows termwise, including `(exp(Z_m)-1)/y` by its positive power series. Term 2 has the stated exact power cancellation. Term 3 requires sqrt(N)>`72.2162666941128`; term 4 requires sqrt(N)>`26.9404922967167` and `35.9206563956223`. Dividing the tail slopes by theta leaves increasing positive denominators. All thresholds hold at N=10000. Independent 201-point logarithmic grids through N=1e12 for c_1/rest and N=1e14 for E_II found zero increases. These grids supplement the analytic arguments.

Evidence: `referee_grid.py/.out`, `referee_local.py/.out`, and the supplied-script rerun outputs `referee_supplied_consts.out`, `referee_supplied_algebra.out`.

## Independent coefficients from the definition

`referee_coefficients.py` evaluates the defining function through its rationalized logit substitution, then extracts coefficients by a Cauchy sum: 2048 angular points, radius 0.96, 65-digit arithmetic. It does not use the proof's M_k contour representation. It returned:

```text
M_0 = -0.553303299720515717370808039043
M_1 =  0.252619096207044794618363520966
```

Let R12 and R13 denote the absolute errors multiplied by N^(3/2)*exp(theta) and N^2*exp(theta), respectively, before division by the theorem constants. Here k=N-1.

| N | M_k | R12 | R13 |
|---:|---:|---:|---:|
| 100 | -1.65220011594776e-10 | 0.594169139411 | 0.232867193260 |
| 107 | -1.70608987364413e-10 | 1.36748240372 | 1.13039344278 |
| 114 | -1.01841660751203e-10 | 1.64712094342 | 1.64770674269 |
| 116 | -8.26583735020128e-11 | 1.62474999943 | 1.69745877265 |
| 158 | 2.43779932331281e-12 | 1.59090765310 | 1.74378702987 |
| 200 | -1.22270461462704e-13 | 1.59721115424 | 1.62331882285 |
| 300 | -2.17221871390157e-16 | 1.12654086420 | 0.938081161232 |
| 381 | 4.58587466651949e-18 | 1.57596334655 | 1.61937319319 |
| 400 | 1.17374036497947e-18 | 1.06547530701 | 1.33489268051 |
| 600 | -3.59803612285491e-22 | 1.26102637861 | 1.11469125087 |
| 1000 | -7.26750314369e-28 | 1.308579915 | 1.18157828 |

All pass both theorems and the main/rest bounds. Among these indices the largest normalized rest is `0.119774723271694153`, at N=107, versus the proved 120000. At N=1000 the rest is approximately `7.1e-38`.

A separate extraction (`referee_convergence.py`) uses 1024 points, radius 0.94, 70 digits, and different x subdivisions. It agrees at N=100,400,600,1000. The largest relative discrepancy is `8.49466668906147e-14`, at N=1000; therefore only about twelve significant digits are treated as checked there. The tiny rest at N=1000 has substantially less relative numerical accuracy than M_k itself. Outputs: `referee_coefficients.out`, `referee_convergence.out`, `referee_convergence_summary.out`.

Additional scaled main-integral checks extend to N=100000, where the complex normalized first-correction error is `1.56489397121694495`, below the analytic bound.

## What was not checked

- I did not reproduce the earlier 85,210-cell certificate for 254.23 or prove the sharper 1030 rest bound. Neither is needed for these two supplied proofs.
- I did not audit the rest of the paper or independently certify its global holomorphy lemma; the coefficient definition and that stated setting were taken as supplied.
- Independent quadratures and coefficient checks are floating-point convergence checks, not interval certificates. I did not numerically check every integer N, or extract M_k from its definition beyond N=1000.
- No fatal defect remains in the audited proof after the explicit local repairs above. This verdict does not validate the uncorrected window inequality or the inaccurate numerical table.
