import mpmath
from mpmath import mp, mpf, sqrt, pi, log, exp, cos, sin, cosh, sinh, nstr
mp.dps = 70
# LSE diamond: 2 cosh(p(x+y)) + 2 cosh(p(x-y)) = e^p  <=> cosh(pX)+cosh(pY) = e^p/2 in X=x+y,Y=x-y
# perimeter = (1/sqrt2) * perimeter of cosh(pX)+cosh(pY) = e^{p(1-ln2/p)}
def P_smsq_level(p, clev):
    # cosh(pX)+cosh(pY) = exp(p*clev), perimeter via arclength in X>=0 quadrant
    E = exp(p*clev)
    yp = lambda X: -sinh(p*X)/sqrt((E-cosh(p*X))**2-1)
    X0 = log(E-1 + sqrt((E-1)**2-1))/p
    return 4*mpmath.quad(lambda X: sqrt(1+yp(X)**2), [0, X0/2, X0])
def P_LSE(p):
    return P_smsq_level(p, 1-log(2)/p)/sqrt(2)
pred = lambda p: 4*sqrt(2) + (4*sqrt(2)*log(2) - 8*log(1+sqrt(2)))/p
for p in [mpf(10), mpf(20), mpf(40)]:
    P = P_LSE(p)
    print(f"p={p}: P={nstr(P,40)}")
    print(f"       pred={nstr(pred(p),40)}  diff={nstr(P-pred(p),6)}")
with open('progress.log','a') as fh:
    fh.write('step_85: softmax square + LSE diamond: single 1/p term, exact constant, NO zeta (parity)\n')
