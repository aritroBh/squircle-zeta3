import sympy as sp
x = sp.symbols('x')
# 1) sech x + tanh x series coefficients delta_d = (-1)^{floor(d/2)} A_d/d!
A = [1,1,1,2,5,16,61,272,1385,7936,50521]
ser = sp.series((1+sp.sinh(x))/sp.cosh(x), x, 0, 11).removeO()
ok1 = all(sp.simplify(sp.expand(ser).coeff(x,d) - (-1)**(d//2)*sp.Rational(A[d], sp.factorial(d)))==0 for d in range(11))
print("sech+tanh coefficients == (-1)^floor(d/2) A_d/d! for d<=10:", ok1)
# 2) d/dx(sech x + tanh x) == (1 - sinh x)/cosh^2 x
ok2 = sp.simplify(sp.diff((1+sp.sinh(x))/sp.cosh(x), x) - (1-sp.sinh(x))/sp.cosh(x)**2) == 0
print("d/dx(sech+tanh) == (1-sinh)/cosh^2:", ok2)
# 3) E(y) duplication identity: Gamma(1+y/2)Gamma(1/2-y/2)/sqrt(pi) == 2^y Gamma(1-y)Gamma(1+y/2)/Gamma(1-y/2)
import mpmath
mpmath.mp.dps = 60
from mpmath import mpf, gamma, sqrt, power
for yy in [mpf('0.23'), mpf('-0.71'), mpf('0.9')]:
    l = gamma(1+yy/2)*gamma(mpf('0.5')-yy/2)/sqrt(mpmath.pi)
    r = power(2,yy)*gamma(1-yy)*gamma(1+yy/2)/gamma(1-yy/2)
    print(yy, mpmath.nstr(l-r, 8))
# 4) eps_m = (-1)^floor((m+1)/2) matches +1,-1,-1,+1
def eps(m): return [1,-1,-1,1][m%4]
print("eps_m == (-1)^floor((m+1)/2):", all(eps(m) == (-1)**((m+1)//2) for m in range(12)))
with open('progress.log','a') as f:
    f.write('step_54: symbolic zigzag identities confirmed; duplication-formula identity E(y) confirmed numerically\n')
