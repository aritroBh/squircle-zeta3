# verify_num: independent high-precision numerics for P4c (Ahmed's Open Problem 4)

Claim checked: `result.md` in this folder,
A(n) = 1 + sqrt2/n + R(n)/n^2 + R3(n)/n^3 + ..., with the stated closed forms for beta, gamma, b3, c3, d3,
and the generalisation E_rho(n) = sqrt(1+rho^2)/(rho n) + R_rho(n)/n^2 + O(log^3 n / n^3).

**Verdict (numerics only):** every asymptotic constant claimed in result.md is confirmed by a free fit to my own
data. Agreement is 47 to 90 digits, depending on the constant. Integer-relation search on the fitted values gives back
exactly the claimed closed forms, including the -zeta(3)/2 in d3, and Euler's gamma does not appear. R_rho is confirmed
to about 75 digits at rho = 1/2 and rho = 3. One factual error sits in the evidence section of result.md: its "55-digit"
E(n) table is correct only to 28-34 significant digits (see section 4). That error does not touch any constant. This
file checks the numbers only. It does not certify that the proof is complete.

All scripts and raw data are in `<scratch>/p4c_num/` (scratch). Nothing else in the project was modified.

## 1. Method (independent of the job's sigma/delta decomposition)

No Phi, no C(n), no hypergeometric function, no Mellin transform, no Tauberian lemma, and none of the job's numbers are
used. I integrate the original integrand directly in the original variables x, y.

Write m = rho*n. The integrand is G = x^(2n-2) y^(2m-2) (n^2 y^2 + m^2 x^2), and I use
sqrt(1+G) - 1 = G/(1+sqrt(1+G)), which has no cancellation.

- **Method 1 (production): polar coordinates about (1,1).** Set x = 1 - r cos(th) and y = 1 - r sin(th), and rescale
  the radius as r = tau/(n k) with k = cos th + rho sin th. Then
  E = int dth (n k)^-2 int tau F dtau.
  - x^(2n-2) and y^(2m-2) are computed as exp((2n-2) log1p(-a)), so the method stays accurate up to n = 1e20.
  - The tau integral uses composite Gauss-Legendre. Unit panels cover |tau - L0| <= 6, where L0 = log(n sqrt(1+rho^2)).
    That is where the boundary layer sits, with branch points at L0 +/- i pi/2. Geometrically growing panels cover the
    rest, and the integral is truncated at tau = L0 + 114, which drops less than e^-228 of the mass.
  - The th integral uses Gauss-Legendre on [0, pi/4] doubled (rho = 1, using the x-y symmetry), or on [0, pi/4] and
    [pi/4, pi/2] (rho not 1).
  - Arithmetic is Arb ball arithmetic (python-flint) at 320 bits for rho = 1 and 360 bits otherwise, with 60x64 or 90x84
    nodes.
- **Method 2 (cross-check): Cartesian.** Set x = 1 - u/n and y = 1 - v/n. The rule is a tensor Gauss-Legendre product
  on unit panels in u and v (panels of width 1/rho in v when rho > 1), at 200 bits with 32 nodes. It is a separate code
  path; the only shared piece is the GL node generator.

Data set: 141 values of n per rho, n = 2*round(10^(3+j/8)/2) for 10^3 <= n <= 10^20, plus 3e3, 3e4, 3e5 and 3e6.
That is 3 x 141 integrals, about 6 CPU-minutes in total.

## 2. Accuracy of my E(n)

Relative change against the next refinement level (Nth, Ntau, bits):

| n | rho | (40,40,200) | (60,64,320) | vs (90,84,360) |
|---|---|---|---|---|
| 1e3 | 1 | -3.6e-58 | -9.0e-92 | ref |
| 1e3 | 1/2 | -2.8e-49 | -4.2e-75 | ref |
| 1e3 | 3 | -3.2e-41 | -1.4e-62 | ref |
| 1e12 | 1 | -4.9e-58 | -1.6e-90 | ref |
| 1e12 | 1/2 | -4.3e-48 | -1.6e-72 | ref |
| 1e12 | 3 | -9.8e-41 | -1.7e-61 | ref |

Production used (60,64,320) for rho = 1 and (90,84,360) for rho = 1/2 and 3, which gives roughly 1e-90 relative.

Method 2 against method 1, as (cart - polar)/polar: -2.7e-56 (n=1e3), -3.6e-56 (1e4) and -2.9e-56 (1e6) for rho = 1;
-5.4e-56 (1e3) and -6.6e-56 (2e3) for rho = 1/2; -3.0e-56 (1e3) for rho = 3. These differences are the 200-bit,
32-node accuracy of the Cartesian run, not disagreement.

Two caveats on this check:
- My first Cartesian run at rho = 3 used unit v-panels. It differed by 7e-34, because the v-singularities sit only
  pi/6 off the real axis. With v-panels of width 1/3 the difference is 3e-56.
- After the fit (section 5), the largest residual is 3e-94 in n*E for rho = 1 and about 1e-96 for rho = 1/2 and 3.
  That is consistent with the data being good to about 90 digits.

## 3. E(n) at the requested n (rho = 1), with the residual chain

r1 = n^2(E - sqrt2/n) - R(n), r2 = n r1 - R3(n), r3 = n r2 - R4(n), where R4 is the formula that result.md §2(d) implies.

| n | E(n) = A(n) - 1 (45 digits) | r1 (~R3/n) | r2 (~R4/n) | r3 (~R5/n) | n r3/(log^5 n/30) |
|---|---|---|---|---|---|
| 1e3 | 0.00138302717521219810699100014938338982746446963 | 0.136027699683 | -0.3445049324 | 0.700772 | 1.337 |
| 3e3 | 0.000466941723005488985728941858089587161630857603 | 0.0675163911622 | -0.1957345196 | 0.443454 | 1.213 |
| 1e4 | 0.000140907818058338745426706214212021723629022561 | 0.0298102605681 | -0.09871466536 | 0.250287 | 1.133 |
| 3e4 | 0.0000470706632683996901344891855535338971992957958 | 0.0136284038477 | -0.05035361951 | 0.140578 | 1.087 |
| 1e5 | 0.0000141344615090732969261614849870854155688137752 | 0.00559298865114 | -0.02304577656 | 0.0709942 | 1.053 |
| 1e6 | 0.00000141410615156815729041413285281391502802510333 | 0.000944145449515 | -0.004669841514 | 0.0170364 | 1.015 |

Further values: E(1e7) = 1.41419922452219489040554717285458489475953623e-7,
E(1e12) = 1.41421356197035338229752588528121991192420151e-12, and
E(1e20) = 1.41421356237309503786747688485234104054045253e-20.

Each residual tends to 0 at the claimed rate. The last column tends to 1, which matches result.md's
R5 ~ (1/30) log^5 n.

## 4. Correction to result.md §4 (precision of the job's E table)

Relative difference (job - mine)/mine for the job's "55-digit" E(n):

| n | 1e3 | 1e4 | 1e5 | 1e6 | 1e7 |
|---|---|---|---|---|---|
| rel. diff | 4.6e-29 | 1.3e-31 | -2.7e-34 | 2.6e-33 | -2.8e-33 |

The job's E(n) is therefore correct to about 28 digits at n=1e3 and 31-34 digits for larger n, not 55. The job's digits
first go wrong at n=1e3: it has ...000149446541..., where the correct value is ...000149383389.... Two independent
quadratures agree on the correct value to 3e-59 at 200 bits.

result.md says these values "match the previous attempt's pkl to all digits", so that pkl carries the same error and is
not an independent check.

Nothing in the asymptotic claims depends on digits beyond about 1e-20, so this is a correction to the evidence, not to
the theorem.

## 5. Free fits (rho = 1)

Model: n E(n) = a0 + sum_{k=1..K} n^-k P_k(log n), with deg P_k = k+1. Every coefficient is free, including sqrt2,
-1/2 and 1/3. The fit is least squares at 170 digits with column scaling.

Fit configurations:

| Config | n range | K | Rows | Parameters | Max residual |
|---|---|---|---|---|---|
| A | 1e3..1e20 | 10 | 141 | 76 | 3e-94 |
| B | 1e4..1e20 | 10 | 132 | 76 | 3e-94 |
| C | 1e5..1e20 | 9 | 123 | 64 | 3e-94 |
| D | 1e6..1e20 | 8 | 114 | 53 | 3e-94 |
| E | 1e3..1e12 | 9 | 77 | 64 | 1e-86 |
| F | 1e3..1e6 only | 4 | 28 | 19 | 1e-28 |

| Coefficient | Fitted (config A, 40 digits) | Closed form in result.md | Diff (A) | Max diff (A-E) | Diff (F: n <= 1e6 only) |
|---|---|---|---|---|---|
| a0 | 1.414213562373095048801688724209698078570 | sqrt2 | 2e-90 | 3e-79 | 4e-23 |
| R: L^2 | -0.5 (exact to 77 digits) | -1/2 | 1e-77 | 2e-65 | 6e-16 |
| beta | -0.653426409720027345291383939270911715962 | lambda/2 - 1 | 9e-76 | 1e-63 | 2e-14 |
| gamma | -2.950163632257010790412668062594455510607 | -(1+sqrt2)+lambda/2-pi^2/12-lambda^2/8 | 2e-74 | 2e-62 | 2e-13 |
| R3: L^3 | 0.3333333333333333333333333333333333333333 | 1/3 | 2e-61 | 5e-52 | 1e-9 |
| b3 | 0.1534264097200273452913839392709117159622 | (1-lambda)/2 | 2e-59 | 5e-50 | 5e-8 |
| c3 | 2.4184737300478041379305747374986031481138 | lambda^2/4-lambda/2+1+zeta2 | 8e-58 | 2e-48 | 8e-7 |
| d3 | 2.4722752293761293983880678485052494425619 | (closed form with -zeta(3)/2) | 1e-56 | 2e-47 | 4e-6 |
| R4: L^4 | -0.125 | -1/8 | 1e-46 | 1e-39 | 5e-5 |
| R4: L^3 | 0.0899534618066529940209746970312108086855 | from §2(d) poly formula | 2e-44 | 2e-37 | 2e-3 |
| R4: L^2 | -1.862142092675839430802239083488496503104 | from §2(d) | 9e-43 | 8e-36 | 3e-2 |
| R4: L^1 | 0.5194479338691649948641392406455201580052 | from §2(d) | 2e-41 | 2e-34 | 0.26 |
| R4: L^0 | -4.973207488596085758846839834986194253850 | from §2(d) | 2e-40 | 1e-33 | 0.84 |
| R5: L^5 | 0.03333333333333333333333333333332978 | 1/30 (result.md, "numerical") | 4e-33 | 3e-28 | 2e-3 |

The zeta(3) term: d3 is pinned to about 1e-47 by every configuration. Dropping the -zeta(3)/2 would shift d3 by 0.601,
so it is excluded by 46 orders of magnitude. With only n <= 1e6 (config F), d3 is still confirmed to 4e-6.

Fitting only the six requested points (1e3, 3e3, 1e4, 3e4, 1e5, 1e6) cannot fix d3. Even with sqrt2, -1/2 and 1/3 fixed,
the five R4 coefficients are unknown. That six-point fit gets beta to 6e-7 and gamma to 8e-6, but b3, c3 and d3 are off
at O(1). This is why I extended the grid.

Structure tests:
- With every P_k raised by one log degree (config A grid, K=9), the extra top coefficients come out as
  L^3 in R: 2e-73, L^4 in R3: 5e-57, L^5 in R4: -2e-42. All are zero, and the lower coefficients are unchanged
  (d3 diff 1e-50).
- An added n^-3/2 term in E gets coefficient 8e-86, which is zero.
- So the expansion has exactly the claimed shape: integer powers of 1/n, and degree k+1 in log n at order n^-(k+1).

As a by-product, the fit gives the R5 L^4 coefficient as -0.0577622650466621091181... = -log(2)/12. It is identified
numerically only; result.md does not claim it.

## 6. Integer-relation identification of the fitted constants

I ran PSLQ at 45 digits on the config-A values. The basis was {1, sqrt2, lambda, lambda^2, lambda^3, zeta2,
zeta2*lambda, zeta3, gamma_E, pi}, so Euler's gamma and pi were both offered as candidates.

- beta  = -1 + lambda/2
- gamma = -1 - sqrt2 + lambda/2 - lambda^2/8 - zeta2/2
- b3    = 1/2 - lambda/2
- c3    = 1 - lambda/2 + lambda^2/4 + zeta2
- d3    = 1 + (3/2)sqrt2 - lambda/2 + lambda^2/8 - lambda^3/24 + zeta2/2 - zeta2*lambda/2 - zeta3/2

These are exactly result.md's closed forms. The coefficients of gamma_E and pi are 0 in every relation, which confirms
the remark that Euler's constant does not appear.

As a separate check, I recomputed R and R3 from result.md's own pieces: R = -sqrt2 - poly1, R3 = 3sqrt2/2 - poly1 + poly2,
with mtilde_k taken from the Taylor coefficients of Gamma((1+z)/2)Gamma(1-z/2)/sqrt(pi). They reproduce the boxed
constants exactly. The R4 implied by §2(d) matches the data to 1e-40 (table above), so the R4 that result.md marks
"numerical" is in fact right to 40 digits.

## 7. Generalisation rho = 1/2 and rho = 3

These are computed directly; I did not use the symmetry E_rho(n) = E_{1/rho}(rho n). For rho = 1/2, n is even, so
m = n/2 is an integer.

The prediction is a0 = sqrt(1+rho^2)/rho and
R_rho = c1/rho - (1/rho)[(L+lr)^2/2 + mt1 (L+lr) + mt2], with lr = log(1+rho^2)/2,
c1 = -(1+rho^3)/(rho sqrt(1+rho^2)), mt1 = 1-lambda and mt2 = 1-lambda+lambda^2/2+zeta2/2.

Free fit, same model as section 5:

| rho | Coefficient | Predicted | Diff A (1e3..1e20, K=10) | Diff B (1e4..) | Diff C (1e5.., K=9) | Diff F (n <= 1e6, K=4) |
|---|---|---|---|---|---|---|
| 1/2 | a0 | sqrt5 = 2.23606797749978969640917366873 | -2e-96 | -5e-96 | -4e-96 | 3e-22 |
| 1/2 | R: L^2 | -1 | 8e-79 | 6e-78 | 3e-78 | -5e-15 |
| 1/2 | R: L^1 | -0.836849190194319136931830847393 | -7e-77 | -5e-76 | -2e-76 | 2e-13 |
| 1/2 | R: L^0 | -6.84493556813006982118981861704 | 1e-75 | 1e-74 | 5e-75 | -2e-12 |
| 3 | a0 | sqrt10/3 = 1.05409255338945977733296451481 | -1e-92 | -1e-92 | -1e-92 | -7e-25 |
| 3 | R: L^2 | -1/6 | 4e-79 | 3e-78 | 1e-78 | 1e-17 |
| 3 | R: L^1 | -0.486048455312359177530587868628 | -4e-77 | -3e-76 | -8e-77 | -5e-16 |
| 3 | R: L^0 | -1.77900671233882915616729397164 | 7e-76 | 5e-75 | 2e-75 | 5e-15 |

Direct residual Rnum - Rpred, with Rnum = n^2(E_rho - a0/n), at the requested n:

| n | E_{1/2}(n) | Rnum - Rpred (1/2) | n*diff/log^3 n | E_3(n) | Rnum - Rpred (3) | n*diff/log^3 n |
|---|---|---|---|---|---|---|
| 1e3 | 0.002176117111772933611692330046930640120886057 | 0.391902 | 1.189 | 0.00104104764149268933531319841738545916225862995 | 0.0444458 | 0.1348 |
| 3e3 | 0.000736750228623598483059510488741964485641113154 | 0.195105 | 1.140 | 0.000349549400518428783540426690235398144350762511 | 0.0210874 | 0.1233 |
| 1e4 | 0.000222613832026122668747739995665065100809282841 | 0.0863988 | 1.106 | 0.000105205404131758250946448584306956689919851156 | 0.00895267 | 0.1146 |
| 3e4 | 0.0000744003693788277407068120095313755176962543333 | 0.0395952 | 1.084 | 0.0000351091983369473622682483108323938700491241125 | 0.00397629 | 0.1089 |
| 1e5 | 0.0000223457787067719214198655960281922094357058239 | 0.0162884 | 1.067 | 0.0000105379790840304991585627854354816129576721243 | 0.00158984 | 0.1042 |
| 1e6 | 0.00000223585870549367601922316717886824168917728452 | 0.00276025 | 1.047 | 0.00000105405224824480153717060319622863921977671035 | 0.000258283 | 0.0979 |

The residuals match result.md §5 (0.392 and 0.086 for rho=1/2; 0.044 for rho=3) and decay like log^3 n / n. The fitted
leading coefficient of R3_rho is 1.000... for rho=1/2, 1/3 for rho=1 and 0.0740740... = 2/27 for rho=3. All three equal
(1+rho)/(6 rho^2). This is my numerical identification; result.md does not claim it. It is a target for anyone writing
out the rho case at the next order.

## 8. What this does and does not establish

**Established numerically, beyond reasonable doubt:**
- the expansion has the claimed shape;
- sqrt2, beta, gamma, b3, c3 and d3 (including -zeta(3)/2) are correct to at least 47 digits;
- R4 from §2(d) is correct to about 33-40 digits, and R5 leads with log^5 n/30;
- R_rho is correct at rho = 1/2 and rho = 3 to about 75 digits.

**Not established here:** that the written proof is complete. Numerics cannot certify these steps:
- the remainder bounds (the Tauberian lemma's O(k!/(3^k N^3)) tail);
- Delta1 = e T1(g) + O(e^2 L^3);
- the termwise Delta2 series.

Also, §5 of result.md gives the rho case only as "same machinery", with no written derivation. If the paper quotes the
rho result as a theorem, that derivation still has to be written out. The numbers support it.

**Correction needed in result.md §4:** its E(n) values are correct to about 28-34 digits, not 55, and the "match" with
the earlier pkl is not independent. The table should be replaced, for example by the values in section 3 above.

## Files (scratch, `<scratch>/p4c_num/`)

- `core.py`: polar integrator (method 1) with Arb arithmetic, plus GL nodes.
- `cart_par.py`, `cart_par_v.py`: Cartesian integrator (method 2); the `_v` version uses v-panels of width 1/rho for rho > 1.
- `conv.py` / `conv.json`: convergence study.
- `prod.py` / `prod.jsonl`: production E values, 3 x 141 n values, 100 digits stored.
- `claims.py`: closed forms, plus R4 and R_rho as result.md's formulas imply them.
- `fit.py`, `fit_all.py`, `fit_rho.py` (+ `.log`, `.json`): free fits.
- `pslq.py`: integer-relation identification.
- `table1.py`: the section-3 table.
