"""Numerical (30 digits) verification of the dual-exponent forms used in the structural discussion.
 P(p)  = perimeter of |x|^p+|y|^p=1                                  (direct arc length)
 (C)   P/8 = int_0^1 (1+u^{p'})^{1/p'} (1+u^2)^{-3/2} du,  p'=p/(p-1)  (Cauchy surface formula: P = int h_K, h_K = ||.||_{p'})
 (F)   P/8 = 1 + eps*M(eps),  M = -F/2, F = int_0^1 [A+B-sqrt(A^2+B^2)] dt,  A=t^{eps-1}, B=(1-t)^{eps-1}, eps=1/p
 (Y)   M(eps) = (1+x) int_0^infty y^{x-1} (1+y^{1+x})^{-1-eps} f(y) dy,  f(y) = (sqrt(1+y^2)-1-y)/2,  x = eps/(1-eps)
 (S)   [analytic, NOT checked numerically] P/8 = (eps/2) * int_0^infty (1+w)^{-1-eps} sqrt(1+w^{2eps-2}) dw, from F = 2/eps - S (two Beta integrals =1/eps each)
precision is set AFTER all imports."""
import mpmath as mp
mp.mp.dps = 30
print("mp.mp.dps =", mp.mp.dps)
def check(p):
    p = mp.mpf(p); eps = 1/p; pp = p/(p-1); x = eps/(1-eps)
    # direct arc length: P = 8 * int_0^{2^{-1/p}} sqrt(1+y'(x)^2) dx  (octant), y=(1-x^p)^{1/p}
    xm = mp.mpf(2)**(-1/p)
    f = lambda s: mp.sqrt(1 + (s**(p-1) * (1-s**p)**(1/p-1))**2)
    Pdir = 8 * mp.quad(f, [0, xm/2, xm])
    Cau = 8 * mp.quad(lambda u: (1+u**pp)**(1/pp) * (1+u*u)**(-mp.mpf(3)/2), [0, 0.5, 1])
    F = mp.quad(lambda t: t**(eps-1) + (1-t)**(eps-1) - mp.sqrt(t**(2*eps-2) + (1-t)**(2*eps-2)), [0, mp.mpf(1)/2, 1])
    Mval = -F/2
    Pf = 8 * (1 + eps*Mval)
    fy = lambda y: (mp.sqrt(1+y*y) - 1 - y)/2
    My = (1+x) * mp.quad(lambda y: y**(x-1) * (1+y**(1+x))**(-1-eps) * fy(y), [0, 1, mp.inf])
    print(f"p={p}: direct {mp.nstr(Pdir,20)} | Cauchy {mp.nstr(Cau*1,20)}*(P/8 form: x8 included) | F-form {mp.nstr(Pf,20)}")
    print("   |direct-Cauchy| =", mp.nstr(abs(Pdir-Cau),3), " |direct-Fform| =", mp.nstr(abs(Pdir-Pf),3), 
          " |M - Y-form| =", mp.nstr(abs(Mval-My),3))
for p in (3, 4, mp.mpf(5)/2):
    check(p)
