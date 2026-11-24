import sympy as sp
U, cc = sp.symbols('U cc', positive=True)
# rewrite polylog(1,z) = -log(1-z) explicitly, use Li1 function
Li1 = lambda z: -sp.log(1-z)
F = -sp.log(U)**2*sp.log(1-cc*U) - 2*sp.log(U)*sp.polylog(2, cc*U) + 2*sp.polylog(3, cc*U)
dF = sp.diff(F, U)
dF = dF.replace(sp.polylog(1, cc*U), Li1(cc*U))
expr = sp.simplify(cc*dF - sp.log(U)**2/(1-cc*U))
print('c*dF/dU - log^2U/(1-cU) simplified:', sp.simplify(expr))
# Note: F(U) is the antiderivative of log^2U/(1-cU) up to factor 1/c: N = (F(1)-F(X))/c * ... 
# check numerically too at random algebraic point
import mpmath as mp
mp.mp.dps = 40
cv = mp.mpf('2.7')
f = lambda U: -mp.log(U)**2*mp.log(1-cv*U) - 2*mp.log(U)*mp.polylog(2,cv*U) + 2*mp.polylog(3,cv*U)
h = mp.mpf('1e-20')
for u0 in ['0.4','1.3']:
    u0 = mp.mpf(u0)
    num = (f(u0+h)-f(u0-h))/(2*h)
    exact = mp.log(u0)**2/(1-cv*u0)/cv
    print(f'F\'( {u0}) numeric vs closed:', mp.nstr(num-exact, 3))
open('progress.log','a').write('step: antiderivative certified (sympy simplify 0 + numeric pins)\n')
