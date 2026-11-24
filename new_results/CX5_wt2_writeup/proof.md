# The weight-two term of every \(M_k\)

**Theorem.** For every \(k\ge2\), the weight-two part of \(M_k\), in the sense of the paper's reduction, is \(T=(\pi^2+6\lambda^2)/24\), where \(\lambda=\log2\).

**Inputs and notation.** We take from the paper \(M=-F/2\), its master integrand
\[
\Phi(\varepsilon,t)=\tfrac12\left(\sqrt{t^{2\varepsilon-2}+(1-t)^{2\varepsilon-2}}
-t^{\varepsilon-1}-(1-t)^{\varepsilon-1}\right),\qquad
M_k=\int_0^1[\varepsilon^k]\Phi\,dt,
\]
and its Euler substitution \(t=(v^2+2v-1)/(4v)\), \(q=\sqrt{t^2+(1-t)^2}=(v^2+1)/(2\sqrt2v)\), \(b=\sqrt2-1\le v\le a=1+\sqrt2\). Put \(\ell_1=\log t\), \(\ell_2=\log(1-t)\), and
\[
D_1=\ell_1'=\frac1{v-b}+\frac1{v+a}-\frac1v,\quad
D_2=\ell_2'=\frac1{v-a}+\frac1{v+b}-\frac1v,\quad
\mathcal F(R)=\frac{R(1/v)}{v^2}.
\]
Hereafter all series are formal in \(\varepsilon\), and \(W^{(k)}=[\varepsilon^k]W\). Canonical keys are \((r,s)\) with \(r\ge s\ge0\), \(d=r+s\ge1\); folding replaces a noncanonical coefficient by \(\mathcal F(R)\) and swaps the indices. Write \(R_{rs}\) for the folded coefficient of \(\ell_1^r\ell_2^s\) in \([\varepsilon^d]\Phi\,dt/dv\). The reducer's triangular system, written explicitly, is
\[
W_{rs}=\varepsilon^dR_{rs}-(r+1)G_{r+1,s}D_1
-[r>s](s+1)G_{r,s+1}D_2-[s=r-1]r\mathcal F(G_{rr}D_1),
\qquad G_{rs}=\mathcal JH(W_{rs}). \tag{*}
\]
Brackets denote indicators. \(H\) retains only principal parts of order at least two at \(\pm i\); \(\mathcal JH\) is their rational primitive normalized to vanish at \(b\). Thus \(G\) is the scripts' `Gm`. This system is just the paper's integration by parts: differentiate \(\ell_1^r\ell_2^s\) and fold. At each fixed level it is uniquely solved downwards in \(d\). Its boundary terms vanish: at \(b\), \(G=O(v-b)\); at \(a\), \(\ell_1^r=O((a-v)^r)\), with \(r\ge1\). These factors dominate any logarithmic divergence. By the reduction's weight convention, simple degree-\(d\) kernels have weight \(d+1\).

The new work below proves the pole bound, the **corrected** leading-coefficient recursion, its initial data and stabilization, and the resulting kernel identity for arbitrary levels. No finite check is used as a proof.

**Lemma A (pole bound).** At level \(k\), for \(d\ge2\), \(H(W_{rs}^{(k)})\) has order at most \(2d-2\) at each of \(\pm i\), and \(G_{rs}^{(k)}\) has order at most \(2d-3\). At degree one, \(H=G=0\).

The initial bound follows from the square-root expansion in Lemma C. The two exponential blocks have no poles at \(\pm i\): for example their first coefficient is \(-(v^2+1)/(2d!v(v-b)(v+a))\). Both \(D_j\) have simple zeros at \(\pm i\), with \(\tau_1=D_1'(i)=(-1-i)/2\), \(\tau_2=D_2'(i)=(-1+i)/2\). Integration of principal parts of order at least two lowers order by one; multiplying by \(D_j\) lowers it again. Normalization adds a constant, and \(\mathcal F\) preserves the paired pole bound. Descending induction in (*) proves the assertion, including analyticity at degree one.

**Lemma B (corrected recursion).** Set \(h_{rs}^{(k)}=[(v-i)^{-(2d-2)}]H(W_{rs}^{(k)})\), \(\rho_{rs}=h_{rs}^{(d)}\), and take absent parents to be zero. For \(k\ge d\ge2\),
\[
h_{rs}^{(k)}=\delta_{kd}\rho_{rs}+
\frac{(r+1)\tau_1h_{r+1,s}^{(k)}
+[r\ge s+2](s+1)\tau_2h_{r,s+1}^{(k)}
+[s=r-1]r\tau_2(h_{rr}^{(k)}-\overline{h_{rr}^{(k)}})}{2d-1}. \tag{B}
\]
A parent top term integrates to \(-h(v-i)^{-(2d-1)}/(2d-1)\); the integration-by-parts minus sign cancels this sign. Lemma A excludes every lower term. At \(i\), inversion has derivative \(1\) and its Jacobian factor is \(-1\), so the folded feed uses minus the conjugate coefficient. At adjacent keys the ordinary second feed and the folded feed are therefore **combined**, not counted twice. This explains the indicator \([r\ge s+2]\).

**Lemma C (initial data).** For \(d=r+s\ge2\),
\[
\rho_{rs}=K_d(-1)^r(1+[r>s])\binom dr,\qquad
K_d=\frac{(-1)^d(2d-3)!!\,i^{d+1}}{2\sqrt2\,d!}. \tag{C}
\]
Indeed, put \(V=(1-t)^2(e^{2x}-1)+t^2(e^{2y}-1)\) and \(U_j=[x^ry^s]V^j\). Expanding the master square root gives its unfolded coefficient
\[
\widetilde R_{rs}=\frac1{2\sqrt2v\,t(1-t)}
\sum_{j=1}^{d}\binom{1/2}{j}q^{2-2j}U_j.
\]
Only \(j=d\) reaches order \(2d-2\). Since each exponential increment starts in degree one, \(U_d=2^d\binom dr(1-t)^{2r}t^{2s}\). At \(v=i\), \(t=(1+i)/2\), so \(U_d=(-1)^ri^d\binom dr\), while \(q'(i)=1/\sqrt2\) and \(t(1-t)=1/2\). Substitution and \(\binom{1/2}{d}=(-1)^{d-1}(2d-3)!!/(2^dd!)\) give \(K_d(-1)^r\binom dr\). Finally \(t(1/v)=1-t(v)\) and \(q(1/v)=q(v)\), so \(\mathcal F(\widetilde R_{sr})=\widetilde R_{rs}\); folding doubles precisely the off-diagonal keys. The exponential blocks add no such poles.

The finite-difference identity used in the earlier derivation is also general:
\[
\sum_{j=0}^d(-1)^j\binom dj j^m(d-j)^{d-m}=(-1)^m d!,\qquad 0\le m\le d.
\]
The polynomial has degree \(d\) and leading coefficient \((-1)^{d-m}\); its signed binomial sum is \((-1)^dd!\) times that coefficient, by the formula for the \(d\)-th forward difference.

**Four identities, with general proofs.** Let \(\widehat h_p=(p+1)\tau_1\rho_{p+1,p}/(4p-1)\). Then
\[
\begin{aligned}
\text{(I1)}\quad &(d+1)K_{d+1}=-i(2d-1)K_d &&(d\ge2),\\
\text{(I2)}\quad &(2d-1)\rho_{rs}=(r+1)\tau_1\rho_{r+1,s}+(s+1)\tau_2\rho_{r,s+1}
&& (r\ge s+2),\\
\text{(I3)}\quad &(4p-3)\rho_{p,p-1}=(p+1)\tau_1\rho_{p+1,p-1}
+p\tau_2(\rho_{pp}-\overline{\rho_{pp}}) &&(p\ge2),\\
\text{(I4)}\quad &\widehat h_p=(1-i)\rho_{pp},\qquad
\widehat h_p-\overline{\widehat h_p}=\rho_{pp}-\overline{\rho_{pp}} &&(p\ge1).
\end{aligned}
\]
(I1) is the ratio of successive expressions in (C). For (I2), both weighted binomials are \((d+1)\binom dr\); use (I1) and \(\tau_2-\tau_1=i\). For (I3), (C) gives \(\overline{\rho_{pp}}=-\rho_{pp}\), \(\rho_{p+1,p-1}=-2p\rho_{pp}/(p+1)\), and \((4p-3)\rho_{p,p-1}=2p\,i\rho_{pp}\). Substitution proves the identity. For (I4), (I1) and \((p+1)\binom{2p+1}{p+1}=(2p+1)\binom{2p}{p}\) give \(\widehat h_p=2i\tau_1\rho_{pp}=(1-i)\rho_{pp}\); its imaginary part is unchanged because \(\rho_{pp}\) is purely imaginary.

**Lemma D (stabilization).** For every \(k\ge r+s\ge2\), off-diagonal keys satisfy \(h_{rs}^{(k)}=\rho_{rs}\). On the diagonal, \(h_{pp}^{(2p)}=\rho_{pp}\) and \(h_{pp}^{(k)}=\widehat h_p\) for \(k\ge2p+1\).

Fix any \(k\) and descend in degree from \(k\), where the data are \(\rho\). Nonadjacent off-diagonal keys follow from (I2). A diagonal key has only its off-diagonal first parent and gives \(\widehat h_p\). An adjacent key uses (I3) if its diagonal parent is initial; otherwise (I4) preserves exactly the difference with its conjugate. These are all cases in (B), proving the assertion for arbitrary \(k\).

**Corollary (the weight-two kernel).** At degree two, Lemma A leaves only double poles and their conjugates. Lemmas C and D give \(h_{20}^{(k)}=-\sqrt2i/4\), \(h_{11}^{(2)}=\sqrt2i/4\), and \(h_{11}^{(k)}=\sqrt2(1+i)/4\) for \(k\ge3\). Integrating \(h/(v-i)^2+\bar h/(v+i)^2\) and subtracting its value at \(b\) therefore gives
\[
\begin{gathered}
G_{20}=\frac{\varepsilon^2}{1-\varepsilon}g_{20},\quad
G_{11}=-\varepsilon^2g_{20}+\frac{\varepsilon^3}{1-\varepsilon}g_{11},\\
g_{20}=\frac{1+\sqrt2}{4}-\frac{\sqrt2}{2(v^2+1)},\quad
g_{11}=-\frac{\sqrt2(v^2+2v-1)}{4(v^2+1)}.
\end{gathered}
\]
The degree-one equation is \(W_{10}=\varepsilon R_{10}-2G_{20}D_1-G_{11}D_2-\mathcal F(G_{11}D_1)\), with the paper's \(R_{10}=(1-1/\sqrt2)/v-2/(v+a)\). Direct substitution gives
\[
(1-\varepsilon)(W_{10}-\varepsilon R_{10})=\varepsilon^2R_{\mathrm{w2}},\qquad
R_{\mathrm{w2}}=\frac1{2v}+\frac1{2(v-a)}-\frac1{v+a}-\frac1{2(v+b)}. \tag{E}
\]
Here is the algebra: \(\mathcal F(D_1)=-D_2\), \(g_{20}+g_{11}=(v^2-2\sqrt2v+1)/(4(v^2+1))\) is invariant under inversion, and \(-2g_{20}D_1+(g_{20}-g_{20}(1/v))D_2=R_{\mathrm{w2}}\), by clearing denominators. These give respectively the zero cubic coefficient and the displayed quadratic coefficient. Thus \([\varepsilon^k]W_{10}=R_{\mathrm{w2}}\) for every \(k\ge2\). It has only simple real poles, so it is exactly the weight-two kernel; no assertion about other weights is needed.

**Evaluation of the common integral.** Independently of the appendix's evaluation, \(R_{\mathrm{w2}}=D_2/2-q'/(1+q)\). Integration by parts uses \(\log((1+q)/2)=0\) at both endpoints and \(\log((1+q)/2)=O(t)\) at \(t=0\), giving
\[
T=\frac{\pi^2}{12}+\int_0^1\frac{\log((1+\sqrt{1-2t(1-t)})/2)}t\,dt
=\frac{\pi^2}{12}-\frac12\operatorname{Li}_2(1/2)
=\frac{\pi^2}{24}+\frac{\lambda^2}{4}.
\]
For the middle equality, differentiation proves \(\log((1+\sqrt{1-z})/2)=-\sum_{n\ge1}\binom{2n}{n}z^n/(2n4^n)\). Insert \(z=2t(1-t)\) and integrate termwise: the beta integral gives \(\binom{2n}{n}\int_0^1t^{n-1}(1-t)^n dt=1/n\). Monotone convergence applied to the negative of the series justifies this step. The last equality is Euler's dilogarithm reflection formula at \(1/2\). This proves the theorem.

**Finite verification only.** `verify_wt2.py` checks the recorded recursion through \(k=12\), exposes the printed recursion's \((k,r,s)=(4,2,1)\) failure, and checks the stated identities over explicit finite ranges. Its rational and binomial/beta certificates also check (E) and the integral evaluation. The unrestricted conclusions above come from the general arguments.
