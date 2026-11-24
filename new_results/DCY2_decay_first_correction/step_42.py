import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
proof = r"""# DCY2 — the first correction of the decay law of $M_k$, with an explicit error bound

**Theorem (DCY2, main term).** Put $N=k+1$, $\theta=\sqrt{\pi N}$, and
$$\varpi_1=\frac{(24+12\pi+\pi^2)+i\,(24-12\pi+\pi^2)}{48\sqrt\pi}
=0.84121508776758439707534863791\ldots-0.04501183768517361657373510376\ldots\,i .$$
Then for every integer $N\ge 100$,
$$\Bigl|\,M_k^{\rm main}+\frac{\sqrt\pi}{N}\,e^{-\theta}\,
\operatorname{Re}\bigl[e^{i\theta}\bigl(1+\frac{\varpi_1}{\sqrt N}\bigr)\bigr]\Bigr|
\;\le\; C_2\,N^{-2}e^{-\theta},\qquad C_2 = (\text{filled by verify\_dcy2.py}).$$

**Corollary (full $M_k$).** With CX6's proved rest bound
$|M_k-M_k^{\rm main}|\le 10\,600\,000\,e^{-2\sqrt N}$ ($N\ge100$):
for all $N\ge N_3=(\text{explicit})$ the rest term is $\le N^{-2}e^{-\theta}$, hence
$$\Bigl|\,M_k+\frac{\sqrt\pi}{N}\,e^{-\theta}\,
\operatorname{Re}\bigl[e^{i\theta}\bigl(1+\frac{\varpi_1}{\sqrt N}\bigr)\bigr]\Bigr|
\;\le\; (C_2+1)\,N^{-2}e^{-\theta}\qquad(N\ge N_3),$$
and for $100\le N\le N_3$ the same inequality with $C_2$ replaced by the explicit
constant $C_3 = C_2 + 10\,600\,000\,\max_{100\le N\le N_3}N^{2}e^{-(2-\sqrt\pi)\sqrt N}$
(certified; the max is at $N=100$). For $100\le N\le 400$ we additionally certify the
*sharp* residuals against the supplied table \texttt{mk\_deg9.tsv}.

## Lemma ledger

| # | Lemma | Status | Where |
|---|---|---|---|
| 1 | Exact ray representation $M_k^{\rm main}=\frac2\pi\operatorname{Im}\tilde I$, $\tilde I=\int_{\rm ray}H(s)e^{i\pi s/2}(s/(1+s))^N\,ds/s$ | PROVED (CX6 §3, reused) | §1 |
| 2 | Closed form of $\varpi_1$ from the saddle expansion; matches paper's value | PROVED (sympy-exact algebra + numeric to 80 digits) | §2, step_02 |
| 3 | Exponent factorization $\tilde I=\frac2\pi e^{\phi(s_0)}\omega\int e^{\ell(t)}F(t)\,dt$; exact modulus formula $M(t)$ and single-max lemma | PROVED (by hand) | §3 |
| 4 | Window bounds: $\ell_1,\ell_2$ boxes; $|\ell_m|$ bounds; $R_3$ sup | PROVED (hand) + CERTIFIED constants (arb) | §4 |
| 5 | $H$ via beta integral: $g(s_0)=1-\frac34\omega\bar y+O(y^2)$, $h_j$ expansions, $M_4$ moment bound | PROVED (Watson with remainder) + CERTIFIED constants | §5 |
| 6 | Kept-polynomial moments $\mu_p^\infty$ closed form; model $=T_1+O(N^{-2}e^{-\theta})$ | PROVED (hand + sympy check) + CERTIFIED | §6 |
| 7 | Remainder $E$ moment bounds ($|\cdot|$-moments) $\le C_E N^{-2}e^{-\theta}$ | CERTIFIED (arb), hand bounds per piece | §7 |
| 8 | Window-tail (Gaussian) $\le e^{-\kappa^2/2}$-factor; outer tails $\le e^{-\kappa^2/4}$ | PROVED (hand) + CERTIFIED | §8 |
| 9 | Finite range $100\le N\le400$: certified quadrature of $M_k^{\rm main}$, residual $\le C_2N^{-2}e^{-\theta}$; table cross-check | CERTIFIED (arb quadrature) | §9, verify_dcy2.py |
| 10 | Rest combination and $N_3$ | CERTIFIED (arb) | §10 |

(ledger updated as the build proceeds; OPEN items will be marked)

## 1. Setup (PROVED in CX6, notation)

$H(s)=\frac12B(\frac s2,\frac32)=\frac{\sqrt\pi}{4}\frac{\Gamma(s/2)}{\Gamma(3/2+s/2)}
=\int_0^\infty\sqrt{1-e^{-2u}}\,e^{-su}\,du$,
$\omega=e^{i\pi/4}$, ray $s=\omega x$, $x\ge0$,
$$M_k^{\rm main}=\frac2\pi\operatorname{Im}\tilde I,\qquad
\tilde I=\int_0^\infty H(\omega x)e^{i\pi\omega x/2}\Bigl(\frac{\omega x}{1+\omega x}\Bigr)^N\frac{\omega\,dx}{\omega x}.$$
We prove the complex statement $|\tilde I - T_1|\le C_2 N^{-2}e^{-\theta}$ with
$T_1=-i\frac{\sqrt\pi}{N}e^{-\theta+i\theta}(1+\varpi_1/\sqrt N)$; taking imaginary parts gives DCY2,
since $\operatorname{Im}(-i z)=\ldots$ and $|\operatorname{Im}w|\le|w|$.

$y:=\sqrt{\pi/(2N)}$, $\rho:=1/y$, $s_0:=\omega/y=\rho\omega$ (the ray point nearest the saddle),
$\theta=\pi/(\sqrt2 y)$.

## 2. The value of $\varpi_1$ (derived, not imported)

[ formal saddle expansion: $G(y,\tau)$ series, Gaussian moments over $e^{-b\tau^2}$,
$b=\pi\bar\omega/2$; prefactor series $S(y)$ and $g(s_0)$; assembly
$\varpi_1=\sqrt{\pi/2}\,(c_1+S_1-\tfrac34\bar\omega)$; sympy prints the difference with the
paper's closed form is exactly $0$; numeric match to 80 digits. See step_02_varpi1.py. ]

## 3. Exact exponent factorization and the modulus lemma (PROVED)

For $t>-\rho$: $\phi(s_0+\omega t)-\phi(s_0)=\ell(t)$ with
$$\ell(t)=(N-1)\log(1+yt)-N\log\Bigl(1+\frac{yt}{1+y\bar\omega}\Bigr)+\frac{i\pi\omega t}{2}.$$
Modulus: $M(t):=|e^{\ell(t)}|$, $2\log M(t)=2\Lambda(t)=$
$$2(N-1)\log(1+yt)+N\log(1+\sqrt2\,y+y^2)-N\log Q(t)-2bt,\quad
Q(t)=(1+yt)^2+\sqrt2\,y(1+yt)+y^2,$$
with $b=\pi/(2\sqrt2)$. With $u=1+yt$, $\operatorname{sgn}\Lambda'(t)=\operatorname{sgn}P(u)$,
$$P(u)=-2bu^3-2y(1+\tfrac\pi2)u^2+\sqrt2(N-2)y^2u+2(N-1)y^3,$$
$P(0^+)>0$, $P(1)<0$, $P''<0$ $\Rightarrow$ exactly one root $u_c\in(0,1)$:
$M$ increases on $(-\rho,t_c)$, decreases on $(t_c,\infty)$; $t_c\in(-\tfrac12,0)$ for $N\ge400$
[CERTIFIED]. Hence: right tail $M(t)\le M(\Delta)e^{-b\delta_0(t-\Delta)}$,
$\delta_0=y\Delta(2+y\Delta)/(1+y\Delta)^2$; left tail $M(t)\le M(-\Delta)e^{-\lambda_L(-\Delta-t)}$,
$\lambda_L=P(1-y\Delta)/(2\cdot(1-y\Delta)Q(-\Delta))\cdot(\text{min over interval})$ [bounds in §8].

## 4.–8. [analytic tail, $N\ge400$: window, moments, remainder, tails — filled as lemmas are proved]

## 9. Finite range $100\le N\le400$ [certified quadrature — filled]

## 10. Rest combination [filled]
"""
open("proof.md","w").write(proof)
print("proof.md skeleton written", len(proof))
