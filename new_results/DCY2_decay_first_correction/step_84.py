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

(Structural form proved; the constant rests on the certified budget evaluated at
grid points and the cell-sweep machinery in `dcy2_sweep.py`; the sweep completes
the all-$N$ certification. Status of every ingredient is in the ledger below.)

**Corollary (full $M_k$).** With CX6's proved rest bound
$|M_k-M_k^{\rm main}|\le 10{,}600{,}000\,e^{-2\sqrt N}$ ($N\ge100$):
$|M_k+\frac{\sqrt\pi}{N}e^{-\theta}\operatorname{Re}[e^{i\theta}(1+\varpi_1/\sqrt N)]|
\le C_3 N^{-2}e^{-\theta}$ for all $N\ge 100$, with $C_3 = C_2 + 10{,}600{,}000\cdot
\max_{100\le N\le N_3} N^2 e^{-(2-\sqrt\pi)\sqrt N}$ (the max is at $N=100$;
$N_3\approx 2.6\times10^4$ is where the rest term drops below $N^{-2}e^{-\theta}$;
sharper values for $N\le 400$ from the certified table `mk_deg9.tsv`).

## Lemma ledger

| # | Lemma | Status |
|---|---|---|
| 1 | Ray representation $M_k^{\rm main}=\frac2\pi\operatorname{Im}\tilde I$, $\tilde I=\int_{\rm ray}H(s)e^{i\pi s/2}(s/(1+s))^N ds/s$, $H(s)=\frac12B(\frac s2,\frac32)$ | **PROVED** (CX6 §3, reused) |
| 2 | $\varpi_1$ closed form from the saddle expansion; equals the paper's 81-digit value | **PROVED** (sympy exact, diff $=0$, two formalisms; numeric 40+ digits) |
| 3 | Exponent factorization $\tilde I=\frac2\pi e^{\phi(s_0)}\omega\int e^{\ell(t)}F(t)dt$; exact modulus $M(t)$; single-max lemma via the cubic $P(u)$ ($P(0^+)>0$, $P(1)<0$, $P''<0$); tail rates $r_\Delta$, $\lambda_L$ | **PROVED** (by hand; rates certified pointwise) |
| 4 | $|\ell_m|\le Ny^{m+1}+y^m/m$; $c_2=-\operatorname{Re}\ell_2\ge0.9b_Ry$; $\operatorname{Re}\ell_1=-y(1+b_Ry/(1+\sqrt2y+y^2))$ (exact) | **PROVED** (hand; $c_2$ is a polynomial inequality certified by interval evaluation) |
| 5 | Watson lemma for $H$: $g(s_0)=1-\frac{3}{4s_0}+O(y^2)$ with $C_g\le1$, $h_j$ bounds, $M_4'$ moment bound; $C_\rho\le 0.5$ | **PROVED** (Watson with remainder) + **CERTIFIED** constants (`step_03_watson.py`) |
| 6 | Kept-polynomial moments $\mu_p^\infty$ in closed form; model $=T_1+O(N^{-2}e^{-\theta})$; the identities $\sum(\text{moment terms})_{O(y)}=0$ and $S_1+g'_1+B_2/(2b)+\gamma_1=\varpi_1\sqrt{2/\pi}$ | **PROVED** (sympy exact, diff $=0$); remainder constants **CERTIFIED** (assembly) |
| 7 | Remainder bounds: $|E(t)|$ by $|\cdot|$-moments with $[T,\Delta]$-split $e^{\gamma|t|^3}$ factors, damped-Gaussian edges (incomplete-$\Gamma$, implementation verified to 12+ digits vs mpmath) | **CERTIFIED** (arb); audit vs truth all OK |
| 8 | Window-tail bounds $\mathrm{TP}_p$ by exact incomplete-$\Gamma$; outer tails via $M(\pm\Delta)$, endpoint monotonicity | **CERTIFIED** (arb) + **PROVED** (monotonicity structure, single-dip $W(u)$) |
| 9 | Budget ratio $E_{\rm total}N^2e^\theta/C_2$ at grid points | **CERTIFIED** (arb, pointwise): $0.033\ (400),\ 0.277\ (2000),\ 0.239\ (10^4),\ 0.074\ (10^5),\ 0.040\ (10^6),\ 0.018\ (10^9),\ 0.015\ (10^{12})$ |
| 10 | Cell sweep $[400,10^5]$ (region B: per-cell factor-assignment); analytic regime $[10^5,\infty)$ (region C: $y$-scaling + log-derivative lemma, $\kappa^2\ge16$) | **CERTIFIED** (`dcy2_sweep.py` / `verify_dcy2.py`) |
| 11 | Finite range $100\le N\le400$: certified quadrature of $M_k^{\rm main}$, residual $\le C_2N^{-2}e^{-\theta}$; cross-check vs `mk_deg9.tsv` | **CERTIFIED** (arb quadrature) |
| 12 | Rest combination and $N_3$ | **CERTIFIED** (arb evaluation) |

## 1. $\varpi_1$ derived, not imported

Two independent saddle-expansion formalisms, both giving the paper's value exactly:

**(a) $\tau$-formalism** (fixed Gaussian $e^{-b\tau^2}$, $b=\pi\bar\omega/2$):
with $G_1=-(1+\pi\bar\omega^2/2)\tau+(\pi\bar\omega/2)\tau^3$,
$G_2=(\tfrac12+\tfrac34\pi\bar\omega^2)\tau^2-(\pi\bar\omega/2)\tau^4$,
$h_1=-3/2$, $h_2=15/8$,
$c_1=\langle G_2+G_1^2/2+h_2\tau^2+G_1h_1\tau\rangle_{e^{-b\tau^2}}
=\frac{(1+i)(15-\pi^2)}{8\sqrt2\,\pi}$ (numeric $0.8312823718875630-0.9364845810788058i$),
and $\varpi_1=\sqrt{\pi/2}\,(c_1+S_1-\tfrac34\bar\omega)$,
$S_1=\pi(1+i)/(6\sqrt2)$. **sympy: diff with the paper's closed form $=0$ exactly.**

**(b) $t$-formalism** (exact $\ell_1,\ell_2$; moment sum $\sum p_p\mu_p$):
the moment-sum $O(y)$ coefficient is **exactly $0$** (sympy), and the prefactor identity
$S_1+g'_1+B_2/(2b)+\gamma_1=\varpi_1\sqrt{2/\pi}$ holds **exactly** (sympy, diff $=0$), with
$B_2=\tfrac12+\tfrac34\pi\bar\omega^2$, $\gamma_1=A_1^2/(4b)$, $A_1=1+\pi\bar\omega^2/2$.

Both were independently **verified numerically**: the $O(y)$ coefficients of the exact
integrals (fitted at $N=10^6$–$10^7$) match the symbolic values to $10+$ digits.

## 2. The analytic tail ($N\ge400$): structure

$y=\sqrt{\pi/(2N)}$, $s_0=\omega/y$, $\ell(t)=(N-1)\log(1+yt)-N\log(1+yt/(1+y\bar\omega))+i\pi\omega t/2$,
$\tilde I=\frac2\pi e^{\phi(s_0)}\omega\int_{-\rho}^{\infty}e^{\ell(t)}F(t)dt$,
$F(t)=H(s_0+\omega t)/H(s_0)$.
Central window $|t|\le\Delta=\kappa\sigma$, $\sigma=(2\cdot0.9b_Ry)^{-1/2}$,
$\kappa^2=\kappa_0^2+4\log(N/400)$, $\kappa_0=2.48$.
Kept polynomial $P(t)$ (14 monomials, $w\le2$ plus $w=4$-odd) with exact shifted-Gaussian
moments $\mu_p^\infty=e^{-\ell_2t_0^2}\sqrt{-\pi/\ell_2}\sum\binom{p}{2m}t_0^{p-2m}\frac{(2m-1)!!}{(-2\ell_2)^m}$.
Remainder $E(t)=F(t)e^{R_3(t)}-P(t)$ bounded by $|\cdot|$-moments (audit ratios $0.001$–$0.7$).
Window-tail by exact incomplete $\Gamma$. Outer tails: $M(t)\le M(\Delta)e^{-r_\Delta(t-\Delta)}$,
$M(t)\le M(-\Delta)e^{-\lambda_L(-\Delta-t)}$ (single-dip lemma: $W(u)=Nu^2(u^2+\sqrt2yu)-(N-1)Q(u)^2$,
$W(1)<0$, $W''>0$ for $u\ge1$ — both hand-proved for all $N\ge100$).

## 3. Certification regions

- **Finite $[100,400]$**: certified arb quadrature of $M_k^{\rm main}$ per integer $N$
  (CX6 `main_arb` template: segments $[0,1,4,16,64,200]$ + proved tail
  $\le\frac{2\sqrt2}{\pi aR^2}e^{-aR}$), check residual $\le C_2N^{-2}e^{-\theta}$;
  cross-check vs `mk_deg9.tsv` (also certified balls).
- **Region B $[400,10^5]$**: per-cell factor-assignment (ratio 1.02 below 3000, 1.3 above):
  decreasing factors at cell-left, increasing at cell-right; products valid.
  Tails by endpoint monotonicity. Each factor's monotonicity proved by hand
  ($y^p$, $S_{\rm sup}$, $g_\pm$, $L_m$, $M_4'$, $\alpha$, $\beta_i$, $|p_p|$,
  $(1-y\Delta)^{-a}$, $B_3$, $\tilde m_j$, $r_\Delta$, $\lambda_L$ at left;
  $\Delta$, $\sigma$, $\kappa$, $\tilde M_j$, $y^{-p}$, $Q_p(\kappa^2/4)$ at right).
- **Region C $[10^5,\infty)$**: every LHS-piece $=y^{a}\cdot(\text{factors})$ with $a\ge3/2$
  after folding $e^{-\kappa^2/4}=\frac{800}{\pi}e^{-\kappa_0^2/4}y^2$; all factors bounded by
  their value at $N_0=10^5$ ($\kappa^2=28.2\ge16$ covers the $Q_p$ polynomial growth;
  $\kappa^3\sqrt y$ decreasing for $N\ge34600$; $\kappa\sqrt y$ decreasing for $N\ge2956$).

## 4. Rest combination

$10{,}600{,}000\,e^{-2\sqrt N}\le N^{-2}e^{-\sqrt{\pi N}}$ iff
$2\log N+\log(10{,}600{,}000)\le(2-\sqrt\pi)\sqrt N$; certified $N_3$
(the function is eventually decreasing; $N_3\approx2.6\times10^4$).
For $100\le N\le N_3$: $C_3 = C_2 + 10{,}600{,}000\cdot N^2 e^{-(2-\sqrt\pi)\sqrt N}$,
max at $N=100$ ($=C_2+1.1\times10^{10}$, crude but explicit);
for $100\le N\le400$ the table gives the sharp residuals.

## 5. Files

`proof.md` (full write-up), `result.md` (this file),
`dcy2_common.py` ($\varpi_1$, model, $\tilde I$ prototype),
`dcy2_budget.py` + `dcy2_audit.py` (bound prototype + audit),
`dcy2_cert.py` (arb bounds, pointwise ratios),
`dcy2_cell.py` (endpoint/region primitives), `dcy2_sweep.py` (cell sweep),
`step_03_watson.py` (certified $C_\rho$), `verify_dcy2.py` (final verifier, <5 min).
"""
open("result.md","w").write(result)
print("result.md written,", len(result), "chars")
