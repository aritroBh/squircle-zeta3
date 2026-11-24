# Check of the Fourier ingredient (iv) used in SQSR's proof: log^2 tan(th) = lam_0 + sum_m lam_m cos(2 m th),
# lam_m = 8 O_m / m (m even), 0 (m odd), O_m = sum_{i<=m/2} 1/(2i-1).  Quadrature, 40 digits.
from mpmath import mp, mpf, quad, log, tan, cos, pi, nstr, linspace
mp.dps = 45
for m in (1, 2, 3, 4, 6, 10, 13, 20):
    v = 4/pi*quad(lambda x: log(tan(x))**2*cos(2*m*x), linspace(0, pi/2, 2*m+2))
    O = sum(mpf(1)/(2*i-1) for i in range(1, m//2+1))
    cl = 8*O/m if m % 2 == 0 else 0
    print(m, nstr(v, 25), ' minus claim:', nstr(v-cl, 3))
