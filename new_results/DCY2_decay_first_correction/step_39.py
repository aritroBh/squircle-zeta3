import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
from mpmath import mp, mpf, pi, sqrt, exp, gamma, beta, log, digamma, psi
mp.dps = 60
from dcy2_common import df2, Hc, ell_coeff, omega, wb, j

N = mpf(400); y = sqrt(pi/(2*N)); s0 = omega/y
H0 = Hc(s0)
# H(s) = sqrt(pi)/4 * Gamma(s/2)/Gamma(3/2+s/2)
# log H = const + logGamma(s/2) - logGamma(3/2+s/2)
# H'/H = (1/2)(psi(s/2)-psi(3/2+s/2));  (H'/H)' = (1/4)(psi1(s/2)-psi1(3/2+s/2)), etc.
def L(m):  # d^m/ds^m log H
    return (mpf(1)/2**(m+1)) * (psi(m, s0/2) - psi(m, mpf(3)/2+s0/2))
# Taylor coeffs h_m = H^{(m)}(s0)/(m! H(s0)) via Faa di Bruno on exp(logH)
import itertools
Ls = [L(m) for m in range(0,6)]
# h_m = sum over partitions of m of prod L_i^{a_i}/a_i! ... use Bell polynomials
from sympy import bell
from sympy.abc import x as sx
h_exact = [mpf(1)]
for m in range(1,5):
    # complete Bell polynomial B_m(L1,...,Lm)/m!
    val = 0
    # expand: h_m = [z^m] exp(sum_{i>=1} L_i z^i / i! * i! ... ) : use recurrence
    # h_m = (1/m) * sum_{k=1..m} L_k * h_{m-k} / (k-1)! ... check: H'/H = L1 ->
    pass
# recurrence: H^{(m)} = sum_{k=0..m-1} binom(m-1,k) L_{k+1} H^{(m-1-k)}  (from H' = H * L1, differentiate)
Hders = [H0]
for m in range(1,6):
    Hders.append(sum(mp.binomial(m-1,k)*Ls[k]*Hders[m-1-k] for k in range(m)))
h_exact = [Hders[m]/(mp.factorial(m)*H0) for m in range(0,6)]
print("h via polygamma:", [mp.nstr(h,12) for h in h_exact])
# fitted t^1 was -0.091927 - 0.001849i ; h1 should be that
print("fitted t1   : -0.091927291861234837 - 0.001848752584128084j")
# check quadrature h1
def HderQ(m):
    return (-1)**m * mp.quad(lambda u: u**m*sqrt(1-exp(-2*u))*exp(-s0*u), [0, mp.inf])
print("quad h1:", HderQ(1)/H0, " quad h2:", HderQ(2)/H0/2)
