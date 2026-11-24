import re
src = open("squircle.py").read()
src = src.replace("tail = m*Jint(mpf(K)/m, p, a) - h(m, K)/2",
                  "tail = power(m,1-s)*Jint(mpf(K)/m, p, a) - h(m, K)/2")
open("squircle.py","w").write(src)
code = r'''
from squircle import *
from mpmath import mp, zeta, dirichlet, mpf, nstr, fabs
setdps(60)
def beta(x): return dirichlet(x, [0,1,0,-1])
def err(a,b):
    d = fabs(a-b)
    return d/fabs(b) if d else mpf('1e-999')
exact = 4*zeta(2)*beta(2)
Z_A = ZA(2, 4, N=30, M=300, K=120, R=12)
print("ZA   :", nstr(Z_A, 50))
print("exact:", nstr(exact, 50))
print("rel err A:", nstr(err(Z_A, exact), 5))
'''
exec(code)
