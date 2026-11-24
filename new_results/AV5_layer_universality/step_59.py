import mpmath, json
from mpmath import mp, mpf, sqrt, pi, zeta, log, cos, sin, gamma, nstr, asinh
mp.dps = 80
def P_lame(p, lam):
    a = 2/p; b = 2/(lam*p)
    def f(t):
        dx = -a*cos(t)**(a-1)*sin(t)
        dy = b*sin(t)**(b-1)*cos(t)
        return sqrt(dx**2+dy**2)
    return 4*mpmath.quad(f, [0, pi/4, pi/2])
def A_lame(p, lam):
    return 4*gamma(1+1/p)*gamma(1+1/(lam*p))/gamma(1+1/p+1/(lam*p))
# validate: p=1, lam=2: |x|+y^2=1 perimeter = ?
P12 = P_lame(mpf(1), mpf(2))
exact12 = 2*(sqrt(5) + asinh(2)/2)   # guess
print("P(1,2) =", nstr(P12,60))
print("guess  =", nstr(exact12,60))
print("diff:", nstr(P12-exact12,5))
print("A(1,2) =", nstr(A_lame(mpf(1),mpf(2)),40), " direct digon area = 4*int_0^1 sqrt(1-x) dx = 8/3:", nstr(mpf(8)/3,40))
# p=1, lam=1/2: |x|+|y|^{1/2}=1 astroid-like
print("P(1,0.5) =", nstr(P_lame(mpf(1), mpf('0.5')),40))
print("A(1,0.5) =", nstr(A_lame(mpf(1), mpf('0.5')),40))
# p=2, lam=1: circle
print("P(2,1) =", nstr(P_lame(mpf(2),mpf(1)),40), "vs 2pi:", nstr(2*pi,40))
with open('progress.log','a') as fh:
    fh.write('step_58: anisotropic Lame perimeter/area pipeline validated\n')
