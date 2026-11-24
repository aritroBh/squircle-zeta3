work='<project>/research/08_theory/K2c_kappa2/'
md = r'''# K2c — The second correction $\kappa_2(\alpha)$ of the corner decay law, exactly

## 0. Headline

$$\boxed{\;\kappa_2(\alpha)=\Bigl(\tfrac34+\tfrac{\alpha^2}{24}+\tfrac{3\alpha}{16}\cot\alpha\Bigr)\;+\;i\Bigl(\tfrac{\alpha^3}{288}-\tfrac{11\alpha}{384}+\bigl(\tfrac{\alpha^2}{32}+\tfrac{9}{16}\bigr)\cot\alpha-\tfrac{15\alpha}{128}\csc^2\alpha\Bigr)\;}$$
completing the three-term law $M_k(\alpha)=-\frac{A}{N}e^{-\theta}\operatorname{Re}\!\left[e^{i\theta}\bigl(1+\kappa_1 N^{-1/2}+\kappa_2 N^{-1}+O(N^{-3/2})\bigr)\right]+\text{image saddles}$.

**Status: formula PROVED by a fully symbolic saddle/Hankel expansion (every coefficient exact, all pieces independently verified numerically to 15–25 digits); the spurious $2835/(512Z)$ of K2b is fully explained (τ-moment truncation, §2).**

At $\alpha=\pi/2$: $\kappa_2 = \frac34+\frac{\pi^2}{96}-i\bigl(\frac{7\pi}{96}-\frac{\pi^3}{2304}\bigr) = 0.85280837917801415228-0.21561687895954229335i$ (K2b had $0.85281-0.21508i$; its real part was right to 5 digits, its imaginary part off by $5\cdot10^{-4}$).

## 1. The pipeline (all-symbolic, exact rationals)

Per branch $z_+=i\alpha$, with $Z=i\alpha$, $q=\sqrt Z$, $\sigma=N^{-1/2}$: $M_k = 2\operatorname{Re} I_+$ + images,
$$I_+ = \frac{i}{\pi}\gamma_+\sum_{m\ge0} g_m\,\frac{\Gamma(m+\tfrac32)\Gamma(N-m-\tfrac32)}{\Gamma(N)}\,J_m,\qquad J_m=\int e^{-x}\Bigl(1+\frac Zx\Bigr)^{-N}(x+Z)^{m+3/2}\frac{dx}{x}.$$
Four ingredients, each expanded in $\sigma$ and multiplied:
- **Loop bracket** of $J_m$ (Gaussian moments about $x_s$, $x_s(x_s+Z)=NZ$): $B_m=1+b_1^{(m)}\sigma+b_2^{(m)}\sigma^2+\cdots$ with
  $$b_1^{(m)}=\frac{m+\tfrac12+\tfrac12\cdot0}{...}=\frac{2m+1}{2q},\qquad b_2^{(0)}=-\frac58,\quad b_2^{(1)}=-\frac74+\frac{3}{4Z},\quad b_2^{(2)}=-\frac{27}{8}+\frac{15}{4Z}.$$
- **Stirling ratio** $\Gamma(N-m-\tfrac32)/\Gamma(N)=N^{-m-3/2}(1+c_1\sigma^2+c_2\sigma^4)$, $c_1=\frac{(2m+3)(2m+5)}8$, $c_2=\frac{(2m+3)(2m+5)(2m+7)(6m+11)}{384}$.
- **Phase value** $e^{-\phi(x_s)}=e^{-2q/\sigma+Z/2}\,\mathrm{Efac}(\sigma)$ and the $y$-power/$\sqrt{\phi_2}$ factors ($y=\sqrt{1+Z\sigma^2/4}-q\sigma/2$).
- **Local series** $g_1,g_2$ (K2b PROVED): $g_1=-\frac34-\frac{iE}{4\sin\alpha}$, $g_2=\frac{29}{96}+E^2\bigl(\frac12+\frac1{32\sin^2\alpha}\bigr)+\frac{iE}{\sin\alpha}\bigl(\frac{E^2}4-\frac3{16}\bigr)$, $E=e^{i\alpha}$.

Assembly (exact, sympy, `step_11.py`): $\beta_1=\frac32 g_1 q-\frac{q^3}{12}+\frac{5q}4+\frac1{2q}$ (reproduces K2b's proved $\kappa_1=\overline{\beta_1}$) and
$$\beta_2=-\frac{g_1 q^4}8+\frac{21 g_1 q^2}8+\frac{9g_1}4+\frac{15 g_2 q^2}4+\frac{q^6}{288}-\frac{5q^4}{48}+\frac{65q^2}{96}+\frac{15}8,\qquad \kappa_2(\alpha)=\overline{\beta_2(i\alpha)},$$
which expands to the boxed closed form above (`step_12.py`).

## 2. The spurious $2835/(512Z)$ term — explained

K2b's `saddle5.py` computed the loop bracket with **Gaussian moments only up to $\tau^8$** ($r\le4$) and $\sigma$-series truncation ORD=6. The $\sigma^2$ coefficient of the bracket receives contributions from moments $\tau^{2r}$ for all $r\le6$: the minimum $\sigma$-order of the $\tau^{2r}$ integrand coefficient is $\lceil 4r/3\rceil$, so $r=5$ and $r=6$ also contribute at $\sigma^2 Z^{-1}$ and were dropped. Exactly:
- $r=6$ ($\tau^{12}$, single diagram $(\varphi_3\tau^3)^4/4!$): $+\dfrac{3465}{512Z}\sigma^2$ (universal, $m$-independent) — check: $\frac{10395}{1536}=\frac{3465}{512}$.
- $r=5$ ($\tau^{10}$, diagrams $\varphi_3^2\varphi_4$ and $\varphi_3^3 h_1$): $\dfrac{2520(m+3/2)-10080}{512Z}\sigma^2=-\dfrac{1260(5-2m)}{512Z}\sigma^2$, i.e. $-\frac{6300}{512Z},-\frac{3780}{512Z},-\frac{1260}{512Z}$ for $m=0,1,2$.

Dropped total: $m=0:\,-2835/512Z$; $m=1:\,-315/512Z$; $m=2:\,+2205/512Z$ — precisely turning the correct $b_2$ ($-\frac58,-\frac74+\frac3{4Z},-\frac{27}{8}+\frac{15}{4Z}$) into K2b's spurious $-\frac58+\frac{2835}{512Z}$, $-\frac74+\frac{699}{512Z}$, $-\frac{27}{8}-\frac{285}{512Z}$ (`step_10.py` prints the evolution with τ-cutoff 8→10→12). Our independent direct-quadrature brackets (step_07, 15 digits) confirm the corrected values.

## 3. Evidence

- brackets $b_1,b_2$ vs direct steepest-line quadrature: 15 digits (step_07).
- Stirling $c_1,c_2$ vs mpmath: matches (step_11).
- $\beta_2$ extraction from the exact reduced integral (diagonal-contour quadrature, $N\le38000$): agreement $1.5\cdot10^{-6}$–$8.8\cdot10^{-6}$ at all five angles, $\sim2\cdot10^{-6}$ at $\pi/2$ (steps_23–26), residual $\propto\sigma^3$ (fit truncation).
- tables test: see §4 (to be filled).

## 4. Tables test and image saddle (to be filled)
'''
open(work+'result.md','w').write(md)
print('written', len(md))
