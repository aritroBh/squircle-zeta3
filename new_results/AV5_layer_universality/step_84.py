import mpmath, json
from mpmath import mp, mpf, sqrt, pi, log, exp, cos, sin, cosh, sinh, nstr
mp.dps = 70
def Fsoft(beta):
    c, s = cos(beta), sin(beta)
    return ((c+s)*log(c+s) + (c-s)*log(c-s))/2 - c*log(c)
JF = mpmath.quad(Fsoft, [0, pi/8, pi/4])
c1_pred = 8*sqrt(2)*JF
print("int_0^pi/4 F d beta =", nstr(JF,50))
print("c1 prediction = 8 sqrt2 J =", nstr(c1_pred,50))
def P_softmax(p):
    y = lambda x: log(exp(p)-cosh(p*x) + sqrt((exp(p)-cosh(p*x))**2-1))/p
    def yp(x):
        v = exp(p)-cosh(p*x)
        return -sinh(p*x)/sqrt(v**2-1)
    x0 = log(exp(p)-1 + sqrt((exp(p)-1)**2-1))/p
    return 4*mpmath.quad(lambda x: sqrt(1+yp(x)**2), [0, x0/2, x0])
for p in [mpf(8), mpf(12), mpf(16), mpf(20)]:
    P = P_softmax(p)
    print(f"p={p}: P={nstr(P,45)}  (P-8)*p={nstr((P-8)*p,30)}")
json.dump({'c1_pred': nstr(c1_pred,60)}, open('ckpt_softmax.json','w'))
