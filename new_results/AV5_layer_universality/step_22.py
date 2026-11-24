import sympy as sp, mpmath, pickle
from mpmath import mp, mpf, sqrt, zeta, log, exp, quad, nstr, factorial
mp.dps = 90
g1 = lambda s: log(1+exp(-2*s))
def Iab(a,b):
    return quad(lambda s: s**a*g1(s)**b, [0,mpf('0.5'),1,2,4,8,16,40])
ck = pickle.load(open('ckpt_integrals90.pkl','rb'))
I  = {eval(k): mpf(v) for k,v in ck['I'].items()}
for key in [(0,8),(2,6)]:
    I[key] = Iab(*key)
    print(key, nstr(I[key],40))
# also verify A7 against direct integration of g7
th = mpmath.tanh
g7 = lambda s: (-mpf(8)/81*s**9*th(s)*(1-th(s)**2) + mpf(4)/9*s**8*(1-th(s)**2)
                + mpf(244)/315*s**7*(th(s)-1))
def eta(n): return (1-mpf(2)**(1-n))*zeta(n)
def A7(a): return -factorial(a+7)*eta(a+8)/mpf(2)**(a+8)*(mpf(8*a)/81*(a+8)+mpf(488)/315)
for a in [0,1]:
    direct = quad(g7 if a==0 else (lambda s: s*g7(s)), [0,1,3,8,20,40])
    print(f"A7({a}): direct {nstr(direct,35)} vs formula {nstr(A7(a),35)}")
pickle.dump({'I':{str(k):nstr(v,80) for k,v in I.items()},
             'K2':ck['K2'],'K3':ck['K3']}, open('ckpt_integrals90.pkl','wb'))
