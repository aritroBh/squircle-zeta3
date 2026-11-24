
from squircle import *
from mpmath import mp, zeta, dirichlet, pi, mpf, nstr, fabs, log10

setdps(60)

def beta(x):
    return dirichlet(x, [0,1,0,-1])

def err(a,b):
    d = fabs(a-b)
    if d == 0: return mpf('1e-999')
    return d/fabs(b)

print("== p=2 anchor, s=4:  Z_2(4) = 4 zeta(2) G ==")
t0 = __import__('time').time()
Z_A = ZA(2, 4, N=30, M=300, K=120, R=12)
t1 = __import__('time').time()
exact = 4*zeta(2)*beta(2)
print("ZA   :", nstr(Z_A, 50))
print("exact:", nstr(exact, 50))
print("rel err A:", nstr(err(Z_A, exact), 5), f" ({t1-t0:.1f}s)")
Z_B = ZB(2, 4, delta=mpf('0.02'))
t2 = __import__('time').time()
print("ZB   :", nstr(Z_B, 50))
print("rel err B:", nstr(err(Z_B, exact), 5), f" ({t2-t1:.1f}s)")
