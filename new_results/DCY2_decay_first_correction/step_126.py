import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
result = r"""# DCY2 — the first correction of the decay law of $M_k$: result

**Theorem (DCY2, main term).** Put $N=k+1$, $\theta=\sqrt{\pi N}$, and
$$\varpi_1=\frac{(24+12\pi+\pi^2)+i\,(24-12\pi+\pi^2)}{48\sqrt\pi}
=0.84121508776758439707534863791-0.04501183768517361657373510376\,i .$$
Then for every integer $N\ge 100$,
$$\Bigl|\,M_k^{\rm main}+\frac{\sqrt\pi}{N}\,e^{-\theta}\,
\operatorname{Re}\bigl[e^{i\theta}\bigl(1+\frac{\varpi_1}{\sqrt N}\bigr)\bigr]\Bigr|
\;\le\; C_2\,N^{-2}e^{-\theta},\qquad C_2 = 1{,}000{,}000 .$$

**All parts certified** by `verify_dcy2.py` (87 s < 5 min):
finite range $[100,400]$ (301 certified quadratures, max res/bound $0.00141$),
region B $[400,10^5]$ (279 cells, worst LHS $=390{,}055$, margin $2.56\times$),
region C $[10^5,\infty)$ (analytic regime, LHS $\le74{,}115$),
part (b) ($C_b=29.3$), rest combination ($N_3=25{,}711$, $C_3=C_2+1.85\times10^{10}$).

**Corollary (full $M_k$).** With CX6's proved rest bound
$|M_k-M_k^{\rm main}|\le 10{,}600{,}000\,e^{-2\sqrt N}$ ($N\ge100$):
for all $N\ge100$,
$$\Bigl|\,M_k+\frac{\sqrt\pi}{N}e^{-\theta}\operatorname{Re}\bigl[e^{i\theta}\bigl(1+\frac{\varpi_1}{\sqrt N}\bigr)\bigr]\Bigr|
\le C_3 N^{-2}e^{-\theta},\qquad C_3=C_2+1.85\times10^{10},$$
and for $N\ge N_3=25{,}711$ the rest term is $\le N^{-2}e^{-\theta}$, so $C_3=C_2+1$ there.
For $100\le N\le400$ the *sharp* residuals are certified against the table `mk_deg9.tsv`
(max res/bound $=0.00141$, i.e. constant $\approx1409$).

## Lemma ledger

| # | Lemma | Status |
|---|---|---|
| 1 | Ray representation $M_k^{\rm main}=\frac2\pi\operatorname{Im}\tilde I$ with $H(s)=\frac12B(\frac s2,\frac32)$ | **PROVED** (CX6 §3) |
| 2 | $\varpi_1$ closed form from the saddle expansion (two formalisms; exact algebra, diff $=0$; matches paper's 81 digits) | **PROVED** (`step_02_varpi1.py`) |
| 3 | Exponent factorization, modulus $M(t)$, single-max lemma ($P(u)$ cubic), tail rates | **PROVED** (by hand) |
| 4 | $|\ell_m|\le Ny^{m+1}+y^m/m$; $c_2\ge0.9b_Ry$; $\operatorname{Re}\ell_1$ exact | **PROVED** (hand; $c_2$ certified) |
| 5 | Watson lemma for $H$: $g(s_0)=1-\frac{3}{4s_0}+O(y^2)$ ($C_g\le1$), $h_j$, $M_4'$; $C_\rho\le0.5$ | **PROVED** + **CERTIFIED** (`step_03_watson.py`) |
| 6 | Kept moments $\mu_p^\infty$ closed form; model $=T_1+O(N^{-2}e^{-\theta})$; identities $\sum(\dots)_{O(y)}=0$ and $S_1+g'_1+B_2/(2b)+\gamma_1=\varpi_1\sqrt{2/\pi}$ | **PROVED** (sympy exact, diff $=0$) + **CERTIFIED** (constants) |
| 7 | Remainder bounds ($|\cdot|$-moments, $[T,\Delta]$-split, damped incomplete-$\Gamma$) | **CERTIFIED** (`dcy2_cert.py`); audit all OK |
| 8 | Window-tail (exact incomplete $\Gamma$); outer tails with the single-dip lemma | **CERTIFIED** + **PROVED** |
| 9 | Region B sweep $[400,10^5]$: 279 cells, worst LHS $=390{,}055$ | **CERTIFIED** (`dcy2_sweep.py`) |
| 10 | Region C $[10^5,\infty)$: LHS $\le74{,}115$ | **PROVED** (y-scaling, §3) + **CERTIFIED** (pointwise) |
| 11 | Finite range $[100,400]$: 301 certified quadratures, max res/bound $0.00141$ | **CERTIFIED** (`step_06_finite.py`) |
| 12 | Part (b): $|{\rm model}-T_1|\le29.3\,N^{-2}e^{-\theta}$; $C_F=10.1$ | **CERTIFIED** (`step_07_partb.py` + assembly) |
| 13 | Rest combination: $N_3=25{,}711$, $C_3$ | **CERTIFIED** (arb) |
| 14 | Validation at $N=100,200,400$ vs `mk_deg9.tsv`: residual $\sim N^{-2}e^{-\theta}$ | **CERTIFIED** (`verify_dcy2.py` §3) |

## 1. $\varpi_1$ derived, not imported

**(a) $\tau$-formalism** (fixed Gaussian $e^{-b\tau^2}$, $b=\pi\bar\omega/2$):
$c_1=\langle G_2+G_1^2/2+h_2\tau^2+G_1h_1\tau\rangle_{e^{-b\tau^2}}
=\frac{(1+i)(15-\pi^2)}{8\sqrt2\,\pi}$
(numeric $0.8312823718875630-0.9364845810788058i$),
$\varpi_1=\sqrt{\pi/2}\,(c_1+S_1-\tfrac34\bar\omega)$, $S_1=\pi(1+i)/(6\sqrt2)$.
**sympy: diff with the paper's closed form $=0$ exactly.**

**(b) $t$-formalism** (exact $\ell_1,\ell_2$; moment sum $\sum p_p\mu_p$): the moment-sum
$O(y)$ coefficient is **exactly $0$** (sympy), and the prefactor identity
$S_1+g'_1+B_2/(2b)+\gamma_1=\varpi_1\sqrt{2/\pi}$ holds **exactly** (sympy, diff $=0$), with
$B_2=\tfrac12+\tfrac34\pi\bar\omega^2$, $\gamma_1=A_1^2/(4b)$, $A_1=1+\pi\bar\omega^2/2$.
Both formalisms were independently **verified numerically** (the $O(y)$ coefficients of the
exact integrals, fitted at $N=10^6$–$10^7$, match the symbolic values to $10+$ digits).

## 2. The analytic tail ($N\ge400$): structure

$y=\sqrt{\pi/(2N)}$, $s_0=\omega/y$, $\ell(t)=(N-1)\log(1+yt)-N\log(1+yt/(1+y\bar\omega))+i\pi\omega t/2$,
$\tilde I=\frac2\pi e^{\phi(s_0)}\omega\int_{-\rho}^{\infty}e^{\ell(t)}F(t)dt$,
$F(t)=H(s_0+\omega t)/H(s_0)$.
Central window $|t|\le\Delta=\kappa\sigma$, $\sigma=(2\cdot0.9b_Ry)^{-1/2}$,
$\kappa^2=\kappa_0^2+4\log(N/400)$, $\kappa_0=2.48$.
Kept polynomial $P(t)$ (14 monomials, $w\le2$ plus $w=4$-odd) with exact shifted-Gaussian
moments $\mu_p^\infty=e^{-\ell_2t_0^2}\sqrt{-\pi/\ell_2}\sum\binom{p}{2m}t_0^{p-2m}\frac{(2m-1)!!}{(-2\ell_2)^m}$.
Remainder $E(t)=F(t)e^{R_3(t)}-P(t)$ bounded by $|\cdot|$-moments with $[T,\Delta]$-split factors
$T=\min(\Delta,c_2/(2\gamma))$ and damped incomplete-$\Gamma$ edges (implementation matches
mpmath to 12+ digits). Outer tails: $M(t)\le M(\Delta)e^{-r_\Delta(t-\Delta)}$,
$M(t)\le M(-\Delta)e^{-\lambda_L(-\Delta-t)}$, single-dip lemma
($W(u)=Nu^2(u^2+\sqrt2yu)-(N-1)Q(u)^2$ convex, $W(1)<0$, $W''>0$ for $u\ge1$).

## 3. Certification regions

- **Finite $[100,400]$**: certified arb quadrature of $M_k^{\rm main}$ per integer $N$
  (segments $[0,1,4,16,64,200]$ + proved tail $\le\frac{2\sqrt2}{\pi aR^2}e^{-aR}$),
  301 values: max res/bound $=0.00141$ at $N=381$.
  Cross-check vs `mk_deg9.tsv` (also certified balls): at $N=100,200,400$,
  $|M_k-M_k^{\rm main}|=1.65\times10^{-13},4.2\times10^{-18},4.7\times10^{-25}$.
- **Region B $[400,10^5]$**: 279 cells (ratio 1.02), each computing the window bracket by
  interval arb and the tails by endpoint monotonicity; the single problematic factor
  $\Lambda(-\Delta)$ enclosed by the mean-value form $\Lambda({\rm mid})+V\cdot{\rm halfwidth}$
  with the algebraically-cancelled derivative $\Lambda'=A_0+N\!A_0'+D_1'$
  ($A_0=\frac12\log(u^2Q(0)/Q)$, $A_0'=[\sqrt2y(u+\sqrt2y)u_N/u-(\sqrt2u+2y)y_N]/(2Q)+(\sqrt2+2y)y_N/(2Q_0)$).
  Every cell LHS $\le C_2=10^6$; worst $=390{,}055$ at $(3075,3137)$.
- **Region C $[10^5,\infty)$**: every LHS-piece $=y^{a}\cdot(\text{factors})$ with $a\ge0$
  after folding $e^{-\kappa^2/4}=\frac{800}{\pi}e^{-\kappa_0^2/4}y^2$; all correction factors
  non-increasing (the $M_4'$, $\alpha^4M_{12}$, $B_6$, $\beta_i$, $E_{\rm drop}$ pieces are $y^0$;
  the $E_{\rm tailmom}$ pieces are $\ge y^0$ with $Q'_p/Q_p<1$ since $Q'_p=Q_p-x^{\deg Q_p}$;
  $\kappa^2=28.2\ge16$ for the $Q_p$ growth; $\kappa\sqrt y$ decreasing for $N\ge635$;
  $\kappa^3\sqrt y$ decreasing for $N\ge1221$; $\lambda_L$ decreasing).
  Pointwise: LHS $=74115,39589,2436,1040,833,803$ at $N=10^5,10^6,10^9,10^{12},10^{16},10^{20}$.

## 4. Part (b) and the model

${\rm model}-T_1=-iye^{-\theta+i\theta}\cdot\delta_F$ with
$|\delta_F|\le C_F=C_S+C_g+C_{S2}+C_E+C_{MS}+\mathrm{cross}\le10.1$:
$C_S\le0.596$, $C_g\le1$, $C_{S2}\le2.543$, $C_E\le3.164$, $C_{MS}\le0.3$, cross $\le2.496$.
$|{\rm model}-T_1|\cdot N^2e^\theta=(\pi^{3/2}/2)(1+\varpi_1/\sqrt N)|\delta_F|\le29.3\,N^{-2}e^{-\theta}$.
Pointwise: $D(N)=1.657,1.620,1.601,1.594,1.578,1.565$ at $N=400,1000,2000,3000,10^4,10^5$.

## 5. Rest combination

$N_3=25{,}711$ where $10{,}600{,}000\,e^{-2\sqrt N}\le N^{-2}e^{-\theta}$ (certified);
for $100\le N\le N_3$: $C_3=C_2+10{,}600{,}000\cdot h(309)$, $h(309)=1748.98$,
so $C_3=1.85\times10^{10}$ (crude but explicit); for $N\le400$ the sharp residuals from
the finite-range quadrature and the table.

## 6. Files

`proof.md` (full write-up), `verify_dcy2.py` (final verifier, 87 s),
`step_02_varpi1.py`, `step_03_watson.py`, `step_06_finite.py`, `step_07_partb.py`,
`dcy2_cell.py`, `dcy2_sweep.py`, `dcy2_regionc.py`,
`dcy2_common.py`, `dcy2_budget.py`, `dcy2_audit.py`, `dcy2_cert.py` (prototype + audit).
All supplied inputs remain read-only.
"""
open("result.md","w").write(result)
print("result.md finalized,", len(result))
