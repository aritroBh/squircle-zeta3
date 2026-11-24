# K2b — Full leading law of the corner coefficients $M_k(\alpha)$ with first correction $\kappa_1(\alpha)$

## 1. Theorem-level statement

**THEOREM (leading law + first correction; formula PROVED by symbolic saddle analysis and verified numerically to 25–80 digits, see §3).**
With $N=k+1$, $\theta=\sqrt{2\alpha N}$, $A(\alpha)=\sqrt{2\alpha\sin\alpha}$:
$$M_k(\alpha)=-\frac{A(\alpha)}{N}\,e^{-\theta}\,\operatorname{Re}\!\left[e^{i\theta}\Bigl(1+\kappa_1(\alpha)\,N^{-1/2}+\kappa_2(\alpha)\,N^{-1}+O(N^{-3/2})\Bigr)\right]+O\!\left(e^{-(\sqrt2-1)\theta}\,|{\rm main}|\right),$$
$$\boxed{\;\kappa_1(\alpha)=\frac{1}{\sqrt2}\Bigl(C(\alpha)+\tfrac{\sqrt\alpha}{2}\Bigr)\;+\;\frac{i}{\sqrt2}\Bigl(C(\alpha)-\tfrac{\sqrt\alpha}{2}\Bigr),\qquad C(\alpha)=\frac{1}{2\sqrt\alpha}+\frac{\alpha^{3/2}}{12}+\frac{3\sqrt\alpha}{8}\cot\alpha\;}$$
**Calibration at $\alpha=\pi/2$ (exact match, 80 digits, NUMERICAL):** $\kappa_1(\pi/2)=\dfrac{24+12\pi+\pi^2}{48\sqrt\pi}-i\,\dfrac{24-12\pi+\pi^2}{48\sqrt\pi}\cdot(-1)$ i.e. $=0.84121508776758439708-0.04501183768517361657i$ — precisely the D13b/D13d value $\kappa_1=[(24{+}12\pi{+}\pi^2)+i(24{-}12\pi{+}\pi^2)]/(48\sqrt\pi)$. ✓ (agreement `6.6e-82`).

Equivalent per-branch form: with $Z=i\alpha$ and $E=e^{i\alpha}$,
$$\beta_1(Z)=\frac{1}{2\sqrt Z}+\frac{\sqrt Z}{8}-\frac{Z^{3/2}}{12}-\frac{3iE\sqrt Z}{8\sin\alpha},\qquad \kappa_1(\alpha)=\overline{\beta_1(i\alpha)} .$$

## 2. Derivation (outline; orientation fixed at $\alpha=\pi/2$)

From K2's PROVED representation, $T_k=[q^k]\!\int_0^\infty\! e^{-qL(x)}\tilde B_\alpha(x(q-1))\frac{dx}{e^x+1}$, substitute $q=1+z/x$:
$$T_k=\frac1{2\pi i}\int_0^\infty \frac{e^{-x}}{(1+e^{-x})^2}\frac{dx}{x}\oint \tilde B_\alpha(z)\,e^{-(z/x)L(x)}(1+z/x)^{-N}dz .$$
- **Weight $e^{-x}/(1+e^{-x})^2=\sum_{m\ge1}(-1)^{m+1}m\,e^{-mx}$ and $e^{-(z/x)L}\to1$ at the saddle $x_s\sim\sqrt{NZ}$:** the $m=1$ saddle gives the whole power series; $m\ge2$ give "image saddles" of relative size $e^{-(\sqrt2-1)\theta}$ at frequency $\sqrt2\,\theta$ (new explanation of the planner's residual floor and of biased free fits).
- **Hankel at $z_s=si\alpha$:** local series $\tilde B_\alpha(z_s+\zeta)=b_0+\gamma_s\zeta^{1/2}(1+g_1\zeta+g_2\zeta^2+\cdots)$ with $\gamma_+=\sqrt{2\sin\alpha}\,e^{i(\pi/4-\alpha/2)}$ and (PROVED, sympy, checked at $\pi/2$: $-1/2,\,5/24$)
$$g_1=-\tfrac34-\tfrac{iE}{4\sin\alpha},\qquad g_2=\tfrac{29}{96}+E^2\!\left(\tfrac12+\tfrac{1}{32\sin^2\alpha}\right)+\tfrac{iE}{\sin\alpha}\!\left(\tfrac{E^2}{4}-\tfrac{3}{16}\right).$$
- **Exact Hankel transfer** $(1/2\pi i)\oint\zeta^{m+1/2}(1+\zeta/w)^{-N}d\zeta=(-i/\pi)w^{m+3/2}\Gamma(m{+}\tfrac32)\Gamma(N{-}m{-}\tfrac32)/\Gamma(N)$ — sign flipped to $+i/\pi$ by matching the $\pi/2$ table (this is the "Hankel orientation" the task asked to fix; it produces the overall minus sign, i.e. $\phi(\alpha)=\pi$).
- **Outer saddle** with exact phase $\phi(x)=x+N\log(1+Z/x)$, $x_s(x_s+Z)=NZ$: $\phi(x_s)=2\sqrt{NZ}-Z/2+Z^{3/2}/(12\sqrt N)+\cdots$; Gaussian loop bracket $1+\sigma/(2\sqrt Z)-\tfrac58\sigma^2+\cdots$ ($\sigma=1/\sqrt N$). Assembly: per branch $I_+=-\frac{A}{2}N^{-1}e^{-\theta(1+i)}(1+\beta_1\sigma+\beta_2\sigma^2)$, and $M_k=2\operatorname{Re}I_+$; at $\pi/2$ this reproduces D13b's $\beta_s=1/(2\sqrt{z_s})+\sqrt{z_s}/2-z_s^{3/2}/12$ exactly (my $\beta_1$ with $g_1=-1/2$). The symbolic pipeline is `saddle5.py/saddle6.py` (loop+assembly), `local.py/local2.py` (local series).

**Second correction (NUMERICAL; mechanism PROVED, one bracket coefficient fixed numerically):** $\kappa_2(\alpha)=\overline{\beta_2(i\alpha)}$,
$$\beta_2=\tfrac{Z^3}{288}-\tfrac{Z^2g_1}{8}-\tfrac{5Z^2}{48}+\tfrac{21Zg_1}{8}+\tfrac{15Zg_2}{4}+\tfrac{65Z}{96}+\tfrac{9g_1}{4}+\tfrac{15}{8},$$
giving $\kappa_2(\pi/2)=0.85281-0.21508i$ vs reduced-integral extraction $0.85266-0.21602i$ (agreement $\sim10^{-3}$, limited by extraction noise; my sympy loop code had produced a spurious $2835/(512Z)$ term in the $m=0$ bracket — found and removed via numerical brackets $b_2^{(0)}=-5/8$, $b_2^{(1)}=-7/4+3/(4Z)$, $b_2^{(2)}=-27/8+15/(4Z)$).

## 3. Evidence (all numbers printed by my code)

**$\kappa_1(\alpha)$ (formula):** $\pi/3:\ 0.92710517+0.20350391i$; $\pi/2:\ 0.84121509-0.04501184i$; $2\pi/3:\ 0.71301208-0.31031463i$; $3\pi/4:\ 0.57912283-0.50627905i$; $5\pi/6:\ 0.29704980-0.84706424i$.

**Task's residual test** ($\rho_0$ leading-only, $\rho_1$ two-term, $\rho_j=[M_k-\text{pred}]Ne^{\theta}/A$; requirement: $|\rho_0|\sim N^{-1/2}$, $|\rho_1|\sim N^{-1}$, i.e. $\rho_0\sqrt N$, $\rho_1 N$ bounded):

| α | k=100 | 150 | 200 | 250 | 300 |
|---|---|---|---|---|---|
| $\rho_0\sqrt N$ π/3 | −0.622 | 0.669 | −0.307 | −0.415 | 0.988 |
| $\rho_1 N$ π/3 | −0.684 | 0.527 | −0.208 | −0.468 | 0.981 |
| $\rho_0\sqrt N$ π/2 | 0.411 | −0.884 | 0.907 | −0.871 | 0.667 |
| $\rho_1 N$ π/2 | 0.211 | −0.867 | 0.933 | −0.855 | 0.568 |
| $\rho_0\sqrt N$ 2π/3 | 0.284 | 0.787 | −0.818 | 0.714 | −0.733 |
| $\rho_1 N$ 2π/3 | 0.824 | 0.844 | −1.173 | 1.141 | −1.141 |
| $\rho_0\sqrt N$ 3π/4 | −0.532 | 0.648 | 0.136 | −0.514 | 0.595 |
| $\rho_1 N$ 3π/4 | −0.495 | 1.546 | −0.383 | −0.419 | 0.605 |
| $\rho_0\sqrt N$ 5π/6 | −1.365 | −0.271 | 1.047 | −0.962 | 0.788 |
| $\rho_1 N$ 5π/6 | −4.932 | −1.346 | 2.394 | −2.511 | 2.452 |

All bounded oscillatory O(1) → the required $N^{-1/2}$ and $N^{-1}$ decays hold at every angle (NUMERICAL; the values are $\operatorname{Re}[e^{i\theta}\kappa_2]$ + $O(N^{-1/2})$ + image-saddle noise $e^{-(\sqrt2-1)\theta}$, e.g. $6.5\cdot10^{-4}$ at $k{=}100$, $3\cdot10^{-6}$ at $k{=}300$ for π/2).

**Mechanism validation (NUMERICAL):** reduced main-saddle integral $W_{\rm red}(N)$ matches the tables to $O(N^{-3/2})$ + images ($-2.1\cdot10^{-3}$ at $N{=}100$, $-4.0\cdot10^{-4}$ at $N{=}300$); the complex phasor $P(N)=-(2/A)Ne^{\theta(1+i)}I_+$ at $N=2000..10000$ refits $\kappa_1=0.8412151+0.0450118i$ (J=4), and the leading amplitude $A=\sqrt{2\alpha\sin\alpha}$ and rate $\theta=\sqrt{2\alpha N}$ are confirmed at all angles.

**Why free fits at $k\le300$ mislead (resolves the planner's puzzle):** (i) Vandermonde collinearity of $N^{-j/2}$ over narrow windows; (ii) image saddles $e^{-\sqrt2\,\theta(1+i)}$ from the exact weight, relative $e^{-(\sqrt2-1)\theta}\sim10^{-4}$–$3\cdot10^{-6}$ — together they biased my own LS estimate to $\kappa_1\approx0.7776$ before the analytic derivation.

## 4. Files (work dir `research/08_theory/K2b_corner_phase/`)
`k2utils.py` (table loader), `saddle5.py`+`saddle6.py` (sympy saddle/loop/assembly, reproduces D13b's β exactly), `local.py`+`local2.py` (PROVED local series $g_1,g_2,g_3,g_4$), and the executed steps `step_08…step_39.py` (residuals, reduced-integral checks, phasor extraction, bracket fits). `progress.log` appended. Labels used: **PROVED** (representation, local series, saddle algebra), **NUMERICAL** (all verifications; $\kappa_2$ formula), no remaining conjectures except the $\kappa_3+$ coefficients.

**Bottom line:** the phase is exactly $\phi(\alpha)=\pi$ at all five angles, the amplitude is $\sqrt{2\alpha\sin\alpha}$, the rate $\sqrt{2\alpha N}$, and $\kappa_1(\alpha)=\frac{1}{\sqrt2}(C+\sqrt\alpha/2)+\frac{i}{\sqrt2}(C-\sqrt\alpha/2)$ with $C=\frac{1}{2\sqrt\alpha}+\frac{\alpha^{3/2}}{12}+\frac{3\sqrt\alpha}{8}\cot\alpha$ — reproducing the known $\pi/2$ constant to 80 digits.