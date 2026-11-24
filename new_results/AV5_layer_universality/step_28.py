import sympy as sp, mpmath, pickle
from mpmath import mp, mpf, sqrt, zeta, log, exp, quad, tanh, nstr, factorial
mp.dps = 90
def eta(n): return (1-mpf(2)**(1-n))*zeta(n)
ck = pickle.load(open('ckpt_integrals90.pkl','rb'))
I  = {eval(k): mpf(v) for k,v in ck['I'].items()}
def K1(m): return 4*factorial(m)/mpf(2)**(m+1)*(eta(m+1)-eta(m))
g1 = lambda s: log(1+exp(-2*s))
g7 = lambda s: (-mpf(8)/81*s**9*tanh(s)*(1-tanh(s)**2) + mpf(4)/9*s**8*(1-tanh(s)**2)
                + mpf(244)/315*s**7*(tanh(s)-1))
# direct ∫ g1 g7
direct = quad(lambda s: g1(s)*g7(s), [0,mpf('0.5'),1,2,4,8,16,32])
# formula
m1,m2,m3 = 9,8,7
Iget = lambda a,b: I[(a,b)]
formula = (mpf(4)/81*(-mpf(m1*(m1-1)*(m1-2))/2*Iget(6,2) + mpf(3*m1)/2*K1(m1-1))
           + mpf(4)/9*(mpf(m2*(m2-1))/2*Iget(6,2) - K1(m2))
           + mpf(244)/315*(-mpf(m3)/2*Iget(6,2)))
print("direct :", nstr(direct,45))
print("formula:", nstr(formula,45))
print("diff   :", nstr(direct-formula,6))
