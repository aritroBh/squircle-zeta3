from mpmath import mp, mpf, mpc, pi, sqrt, exp, cos, sin, gamma, quad
mp.dps = 50
N = 381
sq2 = sqrt(2)
om = (1+1j)/sq2
def integrand(r):
    s = r*om
    return (exp(1j*pi*s/2) * (s/(1+s))**(N-2) / (1+s)**2
            * sqrt(pi)/2 * gamma(s/2+1)/gamma(s/2+mpf(3)/2) * om)
# integrate over [0,200] in segments
pts = [0,1,4,16,64,200]
val = 0
for lo,hi in zip(pts, pts[1:]):
    val += quad(integrand, [lo,hi])
M_main = (2/pi)*val.imag
pi_ = pi
th = sqrt(pi_*N)
varpi_re = (24+12*pi_+pi_**2)/(48*sqrt(pi_))
varpi_im = (24-12*pi_+pi_**2)/(48*sqrt(pi_))
X = (sqrt(pi_)/N)*exp(-th)*( cos(th)*(1+varpi_re/sqrt(N)) - sin(th)*(varpi_im/sqrt(N)) )
print("M_main =", mp.nstr(M_main, 20))
print("X      =", mp.nstr(X, 20))
print("M_main/X =", mp.nstr(M_main/X, 10))
scale = N**2*exp(th)
print("|M - X| N^2 e^th =", mp.nstr(abs(M_main - X)*scale, 10))
print("|M + X| N^2 e^th =", mp.nstr(abs(M_main + X)*scale, 10))
