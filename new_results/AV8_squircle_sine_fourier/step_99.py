import sympy as sp
from mpmath import mp, mpf, pi, log, zeta, pslq, nstr
mp.dps = 45
ln2 = log(2)
m = sp.symbols('m', integer=True, positive=True)
# C_n law (verified n<=17): C_n = (-1)^{(n+1)/2} n (2 floor((n-1)/4)+1) / (2(n^2-1)^2), n odd>=3
# b_n' law: (-1)^{(n-1)/2}/((n^2-1)(2floor((n-1)/4)+1))
# sum over odd n>=3 of b_n' C_n: pair n=4m+1 (m>=1) and n=4m+3 (m>=0)
# n=4m+1: b'=+1/((n^2-1)(2m+1))=1/(8m(2m+1)^2); C=-n(2m+1)/(2(n^2-1)^2)
# n=4m+3: b'=-1/(8(m+1)(2m+1)^2); C=+n(2m+1)/(2(n^2-1)^2)
t1 = (1/(8*m*(2*m+1)**2)) * (-(4*m+1)*(2*m+1)/(2*(8*m*(2*m+1))**2))
t3 = (-1/(8*(m+1)*(2*m+1)**2)) * ((4*m+3)*(2*m+1)/(2*((8*(m+1)*(2*m+1)))**2))
S1 = sp.summation(t1, (m,1,sp.oo))
S3 = sp.summation(t3, (m,0,sp.oo))
SbC = sp.simplify(S1+S3)
print("sum b_n' C_n (n>=3) =", SbC, "=", sp.N(SbC,25))
