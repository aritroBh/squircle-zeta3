import mpmath, json
from mpmath import mp, mpf, sqrt, pi, log, exp, cos, sin, cosh, sinh, tan, nstr
mp.dps = 70
# exact support function on corner sector: h_ca(beta) = sqrt2 cos beta + (sqrt2/p) F(beta)
def Fsoft(beta):
    c, s = cos(beta), sin(beta)
    return ((c+s)*log(c+s) + (c-s)*log(c-s))/2 - c*log(c)
def h_ca(beta, p):
    return sqrt(2)*cos(beta) + (sqrt(2)/p)*Fsoft(beta)
# direct support function by maximization over the curve cosh(px)+cosh(py)=e^p
def y_of_x(x, p):
    z = e := None
    v = exp(p) - cosh(p*x)
    return log(v + sqrt(v**2-1))/p if v > 1 else mpf(0)
def h_softmax(theta, p):
    c, s = cos(theta), sin(theta)
    x0 = log(exp(p)-1 + sqrt((exp(p)-1)**2-1))/p
    def hof(x):
        return x*c + y_of_x(x, p)*s
    # golden section max over [0, x0]
    lo, hi = mpf(0), x0
    gr = (sqrt(5)-1)/2
    x1, x2 = hi-gr*(hi-lo), lo+gr*(hi-lo)
    f1, f2 = hof(x1), hof(x2)
    for _ in range(200):
        if f1 < f2:
            lo, x1, f1 = x1, x2, f2
            x2 = lo+gr*(hi-lo); f2 = hof(x2)
        else:
            hi, x2, f2 = x2, x1, f1
            x1 = hi-gr*(hi-lo); f1 = hof(x1)
    return (f1+f2)/2
p = mpf(10)
print("theta(beta from pi/4), direct h, h_ca formula, diff")
for beta in [mpf('-0.6'), mpf('-0.3'), mpf('0'), mpf('0.3'), mpf('0.6'), mpf('0.7853')]:
    th = pi/4 + beta
    hd = h_softmax(th, p)
    hf = h_ca(beta, p)
    print(f"{nstr(beta,10)} {nstr(hd,30)} {nstr(hf,30)} {nstr(hd-hf,6)}")
