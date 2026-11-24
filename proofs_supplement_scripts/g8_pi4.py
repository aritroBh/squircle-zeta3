"""g8_pi4.py -- the pi^4 part of Theorem 21, exact (sympy).  Input: the two rational functions of n printed by g8_allN.py 4,
  [zeta(4)] P_{4,n} = (n+3)/((2n-1)(2n+1)(2n+3)),   [zeta(2)^2] P_{4,n} = -n/((2n-1)(2n+1)(2n+3)).
Checks: partial fractions; the sums over n of c^{2n}/(2n+s), s=-1,1,3; the closed forms F1, F2; with zeta(4)=pi^4/90 and
zeta(2)^2 = pi^4/36 the combination pi^4/2880 [ (9c - 30/c + 21/c^3) artanh c + 23 - 21/c^2 ] of the statement."""
import sympy as sp
n, c = sp.symbols('n c', positive=True)
f4 = (n + 3)/((2*n - 1)*(2*n + 1)*(2*n + 3)); f22 = -n/((2*n - 1)*(2*n + 1)*(2*n + 3))
print("partial fractions:", sp.apart(f4, n), "  |  ", sp.apart(f22, n))
at = sp.atanh(c)
S = {-1: c*at, 1: at/c - 1, 3: at/c**3 - 1/c**2 - sp.Rational(1, 3)}      # sum_{n>=1} c^(2n)/(2n+s)
F1 = sp.Rational(7, 16)*S[-1] - sp.Rational(5, 8)*S[1] + sp.Rational(3, 16)*S[3]
F2 = -sp.Rational(1, 16)*S[-1] - sp.Rational(1, 8)*S[1] + sp.Rational(3, 16)*S[3]
F1_stmt = (7*c/16 - 5/(8*c) + 3/(16*c**3))*at + sp.Rational(9, 16) - 3/(16*c**2)
F2_stmt = (-c/16 - 1/(8*c) + 3/(16*c**3))*at + sp.Rational(1, 16) - 3/(16*c**2)
print("F1 from partial fractions == closed form:", sp.simplify(F1 - F1_stmt) == 0, "; F2:", sp.simplify(F2 - F2_stmt) == 0)
for s in S:   # the three elementary sums, as power series to c^40
    ser = sp.series(S[s], c, 0, 42).removeO()
    ref = sum(c**(2*k)/sp.Integer(2*k + s) for k in range(1, 21))
    print(f"   sum c^(2n)/(2n+({s})) series to c^40 matches: {sp.expand(ser - ref) == 0}")
for F, f, nm in ((F1_stmt, f4, 'F1'), (F2_stmt, f22, 'F2')):
    ser = sp.series(F, c, 0, 62).removeO()
    ref = sum(f.subs(n, k)*c**(2*k) for k in range(1, 31))
    print(f"   {nm}: Taylor coefficients of the closed form = f(n), n<=30:", sp.expand(ser - ref) == 0)
comb = F1_stmt/90 + F2_stmt/36                       # coefficient of pi^4
stmt = ((9*c - 30/c + 21/c**3)*at + 23 - 21/c**2)/2880
print("F1/90 + F2/36 == (1/2880)[(9c - 30/c + 21/c^3) artanh c + 23 - 21/c^2]:", sp.simplify(comb - stmt) == 0)
