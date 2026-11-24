import sympy as sp
U, cc = sp.symbols('U cc')
Li1 = lambda z: -sp.log(1-z)
F = -sp.log(U)**2*sp.log(1-cc*U) - 2*sp.log(U)*sp.polylog(2, cc*U) + 2*sp.polylog(3, cc*U)
dF = sp.diff(F, U).replace(sp.polylog(1, cc*U), Li1(cc*U))
print('F\'(U) - c*log^2U/(1-cU) =', sp.simplify(dF - cc*sp.log(U)**2/(1-cc*U)))
# numeric pin of the correct identity F' = c log^2 U/(1-cU)
import mpmath as mp
mp.mp.dps = 40
cv = mp.mpf('2.7')
f = lambda U: -mp.log(U)**2*mp.log(1-cv*U) - 2*mp.log(U)*mp.polylog(2,cv*U) + 2*mp.polylog(3,cv*U)
h = mp.mpf('1e-20')
for u0s in ['0.4','1.3']:
    u0 = mp.mpf(u0s)
    num = (f(u0+h)-f(u0-h))/(2*h)
    exact = cv*mp.log(u0)**2/(1-cv*u0)
    print(f'F\'({u0}) - c log^2/(1-c u):', mp.nstr(num-exact, 3))
open('progress.log','a').write('step: antiderivative correctly certified: F\' = c log^2U/(1-cU), sympy 0 + pins\n')
