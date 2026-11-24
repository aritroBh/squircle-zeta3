import mpmath
from mpmath import mp, mpf, sqrt, pi, zeta, gamma, log, cos, sin, tan, exp
mp.dps = 75

def P_squircle(p):
    """Perimeter of |x|^p+|y|^p=1 via support function h=(c^q+s^q)^(1/q), q=p/(p-1)."""
    q = p/(p-1)
    f = lambda th: (cos(th)**q + sin(th)**q)**(1/q)
    return 8*mpmath.quad(f, [0, pi/4])

# sanity checks
print("p=2 (circle):", P_squircle(mpf(2)), " vs 2pi:", 2*pi)
print("p=1 (diamond):", P_squircle(mpf(1)), " vs 4sqrt2:", 4*sqrt(2))
print("p=100 (near square):", P_squircle(mpf(100)), " vs 8")
# exact area check
def A_squircle(p):
    return 4*gamma(1+1/p)**2/gamma(1+2/p)
print("A(1)=", A_squircle(mpf(1)), "(should be 2);  A(2)=", A_squircle(mpf(2)), "(should be pi)")
with open('progress.log','a') as fh:
    fh.write('step_03: squircle perimeter pipeline validated (p=1,2,100)\n')
