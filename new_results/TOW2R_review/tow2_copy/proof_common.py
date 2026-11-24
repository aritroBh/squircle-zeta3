"""Common symbolic setup for the TOW1 identity proofs.
All identities are exact rational-function identities in the symbols (d,p) over Q(i,sqrt2)
(sqrt2 appears only inside explicit constants). A remainder that simplifies to 0 is a
COMPLETE PROOF for all values of the symbols (exact rational arithmetic in sympy)."""
import sympy as sp
i = sp.I
t = { (1,0): (-1-i)/2, (2,0): (-1+i)/2, (1,1): sp.Rational(3,4)*(1-i),
      (2,1): -sp.Rational(3,4)-sp.Rational(3,4)*i, (1,2): 1+sp.Rational(3,4)*i, (2,2): 1-sp.Rational(3,4)*i }
def tc1c(n):
    return (-1)**n * t[(2,n)]   # conj(tau1_n) = (-1)^n tau2_n
def rfun(dd, m, conjg=False):
    s = -2*dd if not conjg else 2*dd
    return 2*(dd-2*m) + s*i - 4*(dd-1)*(dd-2*m)/(2*dd-3)
def R2(dd, m, conjg=False):
    s = -1 if not conjg else 1
    return sp.Rational(1,2)*( (3*dd-3)/2 + s*i*rfun(dd,m,conjg)/2 - (dd-2)**2 )
