# CX6: an explicit rest bound

For every integer **N=k+1 >= 100**,

\[
\boxed{|M_k-M_k^{\rm main}|\le 10\,600\,000\,e^{-2\sqrt N}.}                 \tag{1}
\]

The requested parameters are **C_rest=10,600,000, a=0, c=2>sqrt(pi)**.
Combining (1) with the supplied main-term theorem gives

\[
\boxed{\left|M_k+{\sqrt\pi\over N}e^{-\sqrt{\pi N}}
 \cos\sqrt{\pi N}\right|
 \le1\,210\,000\,000\,N^{-3/2}e^{-\sqrt{\pi N}}.}                         \tag{2}
\]

These constants are deliberately crude. This proves the leading decay law;
it does not certify the first correction varpi_1 in obs:Mk or the sharper
rest exponent sqrt(2*pi).

## Lemma ledger

| Lemma | Status | Location |
|---|---|---|
| Logit form of the definition | **PROVED (by hand)** | §1 |
| Shifted-contour kernel identity, poles and closing arcs | **PROVED (by hand)** | §2 |
| Exchanges and exact closed-m representation | **PROVED (by hand)** | §3 |
| DCY1 endpoint correction and origin residue | **PROVED (by hand)** | §4 |
| Uniform bound including all small-j terms | **PROVED (by hand)** | §§5–6 |
| Main-contour translation segment | **PROVED (by hand)** | §6 |
| Explicit constants and conversion to (2) | **CERTIFIED (arb)** after the hand bounds | §7; verify_rest.py |
| Three table validations, including main-integral error balls | **CERTIFIED (arb)** | §8; verify_rest.py |
| Main-term constant 254.23 for N>=100 | **CERTIFIED (supplied input)** | D13s and D13u; not rerun |
| Sharper rest exponent sqrt(2*pi) and first correction varpi_1 | **OPEN** | Not required for (1) or (2) |

Every new numerical claim is printed by `verify_rest.py`; the complete
output is in `verify_rest.log`. No finite sample substitutes for an all-N
inequality. The proof avoids the infinite m-series entirely.

## 1. Start from the definition

Use M(epsilon)=-F(epsilon)/2 in paper/ieee/main.tex, eq:F and lem:holo,
with the square-root branch specified there. Symmetry restricts its
integral to 0<t<1/2. Substitute g=t/(1-t), then x=-log g. With

\[
 B(z)=1+e^{-z}-\sqrt{1+e^{-2z}},\qquad C(z)=e^zB(z),
\]

ordinary changes of variables in absolutely convergent integrals give,
for real epsilon near zero,

\[
 M(\epsilon)=-\int_0^\infty e^{-\epsilon x}
 (1+e^{-x})^{-1-\epsilon}B((1-\epsilon)x)\,dx
 =-\int_0^\infty e^{-x}(1+e^{-x})^{-1-\epsilon}
 C((1-\epsilon)x)\,dx.                                                    \tag{3}
\]

The differentiation argument of lem:holo identifies [epsilon^k]M with the
paper's M_k. Domination for the new contour formula is proved below.

Put

\[
 H(s)=\int_0^\infty\sqrt{1-e^{-2u}}e^{-su}\,du
 ={1\over2}B(s/2,3/2),\quad K(s)={1\over2\cos(\pi s/2)},\quad
 \omega=e^{i\pi/4},\quad \Gamma:\ s={1\over2}+r\omega\ (0\le r<\infty).
\]

The expression for H follows by substituting exp(-2u) and using
[Euler's beta identity, DLMF 5.12.1](https://dlmf.nist.gov/5.12.E1).
Directly from the integral, for x=Re s>0,

\[
 |H(s)|\le x^{-1},\qquad |H(s)|\le\sqrt{\pi/2}\,x^{-3/2}.                  \tag{4}
\]

The second estimate uses sqrt(1-e^{-2u})<=sqrt(2u) and the Gaussian
integral. H is analytic on Re s>0; its gamma continuation has poles only
at nonpositive even integers.

## 2. Exact shifted-contour identity

For every real z>=0,

\[
 C(z)=1-{2\over\pi}\operatorname{Im}
                   \int_\Gamma H(s)K(s)e^{-sz}\,ds.                       \tag{5}
\]

To prove this, deform Gamma toward the real interval [1/2,2J], indenting
**above** each odd positive integer. Close at Re s=2J with the vertical
segment from 2J to 2J+i(2J-1/2). No pole is inside this contour: the real
cosine poles lie below the indentations and the H poles lie outside.
At s=2n+1 the residue is

\[
 {(-1)^{n+1}\over\pi}H(2n+1)e^{-(2n+1)z}.
\]

Each clockwise upper semicircle contributes -i*pi times this residue.
The straight real pieces have real integrands. On the closing vertical
segment, |K(2J+iy)|=1/(2 cosh(pi*y/2)), |H|<=1/(2J), and |e^{-sz}|<=1.
That integral is O(1/J) and vanishes. The Gamma integral converges
absolutely, as also follows from (13) below. Its imaginary part is thus

\[
 \sum_{n\ge0}(-1)^n H(2n+1)e^{-(2n+1)z}.                                 \tag{6}
\]

This sum is absolutely convergent even at z=0 by (4). The beta recurrence
gives H(1)=pi/4 and H(2n+3)/H(2n+1)=(2n+1)/(2n+4). The binomial theorem
therefore identifies (2/pi) times (6) with
e^z*(sqrt(1+e^{-2z})-1)=1-C(z). This proves (5), including z=0 by dominated
convergence. In particular, all pole contributions and all closing
segments in this derivation have been accounted for.

## 3. Coefficient extraction and the main term

Insert (5) into (3), and define

\[
 J(\epsilon)=\int_0^1(1+t)^{-1-\epsilon}dt
 ={1-2^{-\epsilon}\over\epsilon},\qquad
 \Phi(\epsilon,s)=\int_0^1t^{s(1-\epsilon)}(1+t)^{-1-\epsilon}dt.
\]

Fubini applies absolutely: for |epsilon|<=1/4 and s on Gamma,
Re(s(1-epsilon))>=(1-sqrt(2)/4) Re s>0 and
|(1+t)^{-1-epsilon}|<=1. The inner integrals are uniformly bounded,
and (4), (13) give an integrable outer majorant. The same estimates on a
small complex epsilon disk permit Cauchy coefficient extraction under
both integrals. For real epsilon this yields the identity with Im;
to complexify, write Im as the difference between the upper contour and
its conjugate lower contour first. Equivalently differentiate along the
real epsilon axis, justified by the same uniform bounds.

Set

\[
 J_k={(-1)^k(\log2)^{k+1}\over(k+1)!},\qquad
 P_k(s)={(-1)^k\over k!}\int_0^1{t^s\over1+t}
                        (s\log t+\log(1+t))^kdt.
\]

The exact, absolutely convergent representation is

\[
 \boxed{M_k=-J_k+{2\over\pi}\operatorname{Im}
                         \int_\Gamma H(s)K(s)P_k(s)\,ds.}                \tag{7}
\]

Deleting (1+t)^{-1-epsilon} in Phi gives 1/(1+s(1-epsilon)), whose kth
coefficient is Q_k(s)=s^k/(1+s)^{k+1}. This is the m=0 contribution.
Put A_k=P_k-Q_k and E(s)=exp(i*pi*s/2), the n=0 contribution to K.

The question's M_main is (2/pi)Im integral_ray H E Q_k ds, since
Q_k=(s/(1+s))^N/s. Rotation from the positive real axis to that ray is
valid by Cauchy's theorem: H E Q_k is analytic in this sector, its origin
singularity is removable for k>=1, and the pole at -1 lies outside.
At zero it is O(s^{k-1}); on a large arc it is O(R^{-5/2}) by (4), since
|E|<=1 and |s/(1+s)|<=1. Thus both the small and large arcs vanish.

Translation from the ray to Gamma crosses no poles and retains the real
segment [0,1/2]. Its imaginary contribution is **not** zero for H E Q_k.
Consequently the exact rest decomposition is

\[
 \begin{split}
 M_k-M_k^{\rm main}={}&-J_k+{2\over\pi}\operatorname{Im}
 \left\{\int_\Gamma HKA_k\,ds+
 \int_\Gamma H(K-E)Q_k\,ds-\int_0^{1/2}HEQ_k\,ds\right\}.
 \end{split}                                                             \tag{8}
\]

The connecting segments at infinity vanish by (4) and exponential decay.

## 4. What DCY1 step_50 misses

H(s)=1/s+O(1), K(0)=1/2, and P_k(0)=J_k. The origin residue of H K P_k
is therefore J_k/2. Between Gamma and the unshifted ray there are no poles;
the strip meets the real axis only between 0 and 1/2. H K P_k is real
there. The small arc of angle pi/4 gives

\[
 \operatorname{Im}\int_{\rm ray}HKP_k ds
 =\operatorname{Im}\int_\Gamma HKP_k ds-\pi J_k/8.
\]

Only the imaginary part of the ray integral is intended: its real part
has a logarithmic divergence at zero. Formula (7) becomes

\[
 M_k=-\tfrac34 J_k+{2\over\pi}\operatorname{Im}
                             \int_{\rm ray}HKP_k ds.                     \tag{9}
\]

Step_50 uses -3/4*(-log 2)^k/k! instead of -3/4*J_k. The missing factor
is **log(2)/(k+1)**. Its separate n-integrals also require care. Since
P_k(0)=J_k is nonzero, their scaled imaginary parts tend to J_k/2, rather
than zero, as n tends to infinity. To see this subtract J_k*e^{-s}/s
from H(s)P_k(s). The remainder is bounded at zero and its integral against
exp(i*(2n+1)*pi*s/2) tends to zero. The subtracted term has imaginary
integral tending to pi*J_k/4 on the pi/4 ray.

Thus that alternating n-series is **not ordinarily convergent**. It can
be read in the Abel sense (replace (-1)^n by (-tau)^n and let tau increase
to one), or avoided by using (7). For the Abel limit, the summed factor
is E/(1+tau*E^2); its imaginary integrand has a uniform bound near zero
after subtracting the real 1/s part, and is uniformly integrable at
infinity. Formula (9) uses the summed cosine kernel, without truncating
this series. The verifier prints both elementary terms at k=8,10,12.

In the Laguerre route, the coefficient extraction before using cut
formulas also includes j=0: the inner index runs to min(m,k), not just
min(m,k-1). We do not assume the j>=1-only display in D13d is exact
without its endpoint correction. Formula (7) accounts for the endpoint
directly. Also the correct coefficients are
c_i^(m)=(-1)^m e_i(1,1/2,...,1/m); several exploratory DCY1 scripts use
the reversed index e_{m-i}. Neither infinite m-summation is used here.

## 5. Bound the closed kernel uniformly

Write beta(L)=log(1+exp(-L)), with its analytic branch on Re L>0, and
rho=|s|. Substituting L=-log t gives

\[
 A_k(s)={1\over k!}\int_0^\infty e^{-(s+1)L}
 \left({(sL-\beta(L))^k\over1+e^{-L}}-(sL)^k\right)dL.                   \tag{10}
\]

For |arg L|<=pi/4 and 0<=v<=1,

\[
 |1+v e^{-L}|\ge\tfrac12,\qquad |\beta(L)|\le2e^{-\operatorname{Re}L}.
                                                                            \tag{11}
\]

Indeed L=x+iy with |y|<=x. For x<=log 2, cos y>=0 and the first modulus
is at least one. For x>=log 2, the triangle inequality gives 1/2.
Integrate e^{-L}/(1+v e^{-L}) in v to obtain the second bound. We use
1/2<=log 2<1<pi/2, from the integral definition of log and pi>3.

Expand the kth power in (10), a **finite** sum. For each i>=1 term rotate
the L contour clockwise to arg L=-arg(s+i+1); for the denominator
correction (i=0) rotate to -arg(s+2). All these angles lie in [-pi/4,0].
There are no zeros of 1+exp(-L) in the rotating wedge, and its logarithm
has the branch fixed above, so no pole or branch cut is crossed. Small
arcs vanish because the integrands are bounded times nonnegative powers
of L. Large arcs vanish by (11): after extracting exp(-i L) from beta^i,
the decay is controlled by Re((s+1+i)L)>0 throughout the rotation.
The denominator correction similarly extracts exp(-L).

On the final ray the exponential rate is |s+1+i|, or |s+2|. The monomial
integral is a factorial and cancels its finite-binomial coefficient,
giving the bound

\[
 \boxed{|A_k(s)|\le {2\rho^k\over|s+2|^{k+1}}
       +2\sum_{i=1}^k{2^i\over i!}
                    {\rho^{k-i}\over|s+i+1|^{k-i+1}}.}                   \tag{12}
\]

This includes i=k and has no infinite m-sum. The factorial i! is the
smallness mechanism for small j=k-i.

Let a=pi/(2 sqrt(2)). On Gamma,

\[
 |H(s)|\le2,\qquad |K(s)|\le2e^{-ar},\qquad |K(s)-E(s)|\le2e^{-3ar}.
                                                                            \tag{13}
\]

For the latter two inequalities put y=2ar. The denominator
1+exp(i*pi*s)=1+i exp((-1+i)y) has modulus at least 1/2: for y>=log 2
use the triangle inequality, and for y<=log 2 its squared modulus is
(exp(-y)-sin y)^2+cos^2 y >=cos^2(log 2)>1/4. Here log 2<1<pi/3.
Now use K=E/(1+E^2) and K-E=-E^3/(1+E^2).

For integers j>=0, real q>=1, and lambda>0, put

\[
 F(j,q;\lambda)=\int_0^\infty e^{-\lambda r}
                     {\rho^j\over|s+q|^{j+1}}dr.
\]

Since 0<=arg s<=pi/4, |s+q|>=rho+q/sqrt(2), rho<=r+1/2, and |s+q|>=1.
Use log(1+v)>=v/(1+v), and set t=r+1/2+q/sqrt(2). Reserve
exp(-lambda*r/4) and apply AM–GM to (3 lambda/4)t+jq/(sqrt(2)t).
Integration of the reserved exponential gives

\[
 F(j,q;\lambda)\le {4\over\lambda}
 \exp\left({\lambda\over2}+{\lambda q\over\sqrt2}\right)
 \exp\left(-\sqrt{{3\lambda jq\over\sqrt2}}\right).                       \tag{14}
\]

The translation actually gives the smaller positive prefactor
exp((3 lambda/4)*(1/2+q/sqrt(2))); we enlarged it in (14).

## 6. Assemble the all-N estimate

Define

\[
 x=2e^{\pi/4},\qquad K_0={64\over\pi a}e^{a/2+\pi/4},\qquad
 T_0={32\over3\pi a}e^{3a/2+3\pi/4}.
\]

Equations (12)–(14) imply

\[
 {2\over\pi}\left|\int_\Gamma HKA_k ds\right|
 \le K_0\left[e^{\pi/4}e^{-\sqrt{(3\pi/2)k}}
 +\sum_{i=1}^k{x^i\over i!}
           e^{-\sqrt{(3\pi/4)(i+1)(k-i)}}\right].                         \tag{15}
\]

For 1<=i<=floor(k/2),
(i+1)(k-i)-2(k-1)=(i-1)(k-i-2)>=0. Thus this product is at least
2(N-2). Since pi>3 and N>=100, (3*pi/2)*(N-2)>=4N. These exponential
factors, and the separate i=0 factor in (15), are at most exp(-2 sqrt(N));
their coefficient sum is at most exp(x).

For i>k/2 discard the saddle exponential and bound the exponential series:

\[
 \sum_{i>k/2}{x^i\over i!}
 \le e^{-k/2}\sum_{i\ge0}{(ex)^i\over i!}
 \le e^{ex}e^{-2\sqrt N}.                                                 \tag{16}
\]

Here k/2=(N-1)/2>=2 sqrt(N) for N>=100; for example
(N-1)/2>=N/4>=2 sqrt(N). Equation (16) handles **all small-j terms**,
including j=0, without a fixed-m uniformity assumption.

The non-main cuts of the m=0 kernel satisfy

\[
 {2\over\pi}\left|\int_\Gamma H(K-E)Q_k ds\right|
 \le T_0e^{-\sqrt{(9\pi/4)k}}\le T_0e^{-2\sqrt N}.                        \tag{17}
\]

This follows from (14) with lambda=3a, q=1, j=k; the final comparison
uses pi>3 and N>=100.

Two elementary contributions remain in (8). First,
|J_k|=(log 2)^N/N!<=2^{-(N-1)}<=exp(-(N-1)/2)<=exp(-2 sqrt(N)).
Here N!>=2^{N-1} and log 2>=1/2 follow by hand. Second, on the real
translation segment use H(s)<=1/s and substitute u=s/(1+s):

\[
 {2\over\pi}\left|\int_0^{1/2}HEQ_k ds\right|
 \le {2\over\pi}\int_0^{1/2}{s^{k-1}\over(1+s)^{k+1}}ds
 ={2\over\pi k}3^{-k}\le e^{-2\sqrt N}.                                  \tag{18}
\]

For the last comparison, 2/(pi*k)<=1 and
3^{-k}<=exp(-k/2)<=exp(-2 sqrt(N)). Therefore (8) gives, for all N>=100,

\[
 |M_k-M_k^{\rm main}|\le C_*e^{-2\sqrt N},\qquad
 C_*=2+T_0+K_0(e^{\pi/4}+e^x+e^{ex}).                                    \tag{19}
\]

All range reductions in this section are hand inequalities, with no
numerical subdivision in N. Bounds pi>3 and pi<4 follow from elementary
inscribed/circumscribed geometry; Arb also prints pi and sqrt(pi).

## 7. Certified constants and conversion

`verify_rest.py` evaluates (19) with outward rounding using python-flint
arb at 90 decimal digits. The full enclosures, including positive margins
to the selected integer constants, are in verify_rest.log. In particular
C_*<10,600,000, proving (1). Printed enclosures include:

| Quantity | Arb enclosure (rounded outward for display) |
|---|---|
| K_0 | [70.0985853351116237515327912023 +/- 1.45e-29] |
| T_0 | [170.657769858025927890213093913 +/- 2.52e-28] |
| C_* | [10578907.8235293971071051548801 +/- 4.83e-23] |
| 2-sqrt(pi) | [0.227546149094483972701832516659 +/- 1.46e-31] |

For delta=2-sqrt(pi)>0, the hand maximization

\[
 \sup_{u>0}u^3e^{-\delta u}=(3/(e\delta))^3
\]

implies that the additional coefficient in (2) is at most
10,600,000*(3/(e*delta))^3. The verifier certifies
254.23+10,600,000*(3/(e*delta))^3<1,210,000,000, proving (2) conditional
only on the main-term theorem supplied in the question. Its printed
enclosures are [114.096375458588480449897257230 +/- 1.64e-28] for the
conversion factor and [1209421834.09103789276891092664 +/- 3.74e-21]
for the resulting unrounded total constant.

## 8. Validation and reproduction

The supplied mk_deg9.tsv has **no textual header**. Its rows are an index,
an M_k ball and a second very small ball. The verifier parses the complete
second column with arb, preserving its supplied radius; the third column
is not silently reinterpreted as an error bound.

M_main is evaluated on the unshifted ray, which equals the specified
positive-real integral by §3. Mpmath uses extra working precision and
prints 60 significant digits. Independent Arb integration on 0<=r<=200
encloses the same main integral, after cancelling its removable pole:

\[
 e^{i\pi s/2}\left({s\over1+s}\right)^{N-2}{1\over(1+s)^2}
 {\sqrt\pi\over2}{\Gamma(1+s/2)\over\Gamma(3/2+s/2)}\,\omega.
\]

Only integer powers and meromorphic functions are used, so the Arb
analytic callback has no unhandled log or square-root branch. The tail
bound follows from |H(r*omega)|<=sqrt(2)/r and |s/(1+s)|<=1:

\[
 {2\over\pi}\left|\int_R^\infty HEQ_k\omega\,dr\right|
 \le {2\sqrt2\over\pi aR^2}e^{-aR},\qquad R=200.                          \tag{20}
\]

The verifier prints this allowance and adds its upper endpoint as a
symmetric error radius. It then encloses abs(table M_k - M_main) and
asserts this is below (1) at k=100,200,400. The validation is certified
by the Arb balls, not by agreement between two floating-point answers.

| k | Certified abs(M_k-M_main), outward display enclosure | Bound (1), outward display enclosure |
|---|---|---|
| 100 | [1.69252447066814635e-13 +/- 1.77e-31] | [0.0197740128995278421 +/- 2.65e-21] |
| 200 | [3.81140793878010622e-18 +/- 4.66e-36] | [5.13963137274728551e-6 +/- 1.08e-24] |
| 400 | [4.79912536555363781e-25 +/- 3.73e-43] | [4.28376284616536664e-11 +/- 4.65e-29] |

The log retains the much narrower full integration balls. A separate
check of the new exact representation (7), at k=0, uses
P_0(s)=(psi(1+s/2)-psi((1+s)/2))/2. This follows by expanding 1/(1+t),
pairing even and odd terms, and applying the defining digamma difference
series. Its certified result minus the table is [+/- 9.36e-81], containing
zero. The added tail allowance is 8*exp(-200*a)/(pi*a), from (13) and
|P_0(s)|<=1. This checks the elementary correction and contour orientation
independently of the very loose rest inequality.

The complete recorded run passed in 13.356 seconds.

Run from this directory:

```sh
'<project>/.venv/bin/python' verify_rest.py
```

All supplied inputs remain read-only. The existing main-term certificate
is assumed as requested; its cell sweep is not repeated.
