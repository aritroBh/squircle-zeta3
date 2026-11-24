# step_47: pin down beta for the golden Li3 classes
import mpmath as mp
mp.mp.dps = 150
def Li3(x): return mp.re(mp.polylog(3, mp.mpf(x)))
def Li2(x): return mp.re(mp.polylog(2, mp.mpf(x)))
def L3(x):
    x = mp.mpf(x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(1-x)
def chi5(n):
    r = n % 5
    return 0 if r==0 else (1 if r in (1,4) else -1)
L5 = mp.nsum(lambda n: chi5(int(n))/mp.mpf(n)**3, [1, mp.inf])
z3 = mp.zeta(3)
phi = (1+mp.sqrt(5))/2
X1, Xm1 = L3(1/phi), L3(-1/phi)
alpha = (X1+Xm1)/(2*z3)
beta  = (X1-Xm1)/(2*L5)
import sympy as sp
print("alpha:", sp.nsimplify(str(alpha)))
print("beta  =", mp.nstr(beta, 50))
# try beta as rational: use pslq [beta, 1]
rel = mp.pslq([beta, mp.mpf(1)], tol=mp.mpf('1e-120'), maxcoeff=10**20)
print("pslq [beta,1]:", rel)
# also check X1 - Xm1 against q*L5 with rational q:
q = (X1-Xm1)/L5
rel2 = mp.pslq([q, mp.mpf(1)], tol=mp.mpf('1e-120'), maxcoeff=10**20)
print("pslq [(X1-Xm1)/L5, 1]:", rel2)
# check: maybe {1/phi} = (1/10){1} + beta*{chi}: numeric value of beta*L5*2:
print("2*beta*L5 =", 2*beta*L5)
print("X1-Xm1    =", X1-Xm1)
