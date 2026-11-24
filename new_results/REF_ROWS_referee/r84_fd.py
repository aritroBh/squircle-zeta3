# Finite-difference check of P'''(1) for the octagon directly from P(p) (no termwise differentiation).
from mpmath import mp, mpf, sqrt, log, cos, sin, pi, quad, nstr
mp.dps = 40
m = 4
def P(p):
    def f(th):
        S = Sp = mpf(0)
        for i in range(1, m+1):
            c = cos(th-i*pi/m); s = sin(th-i*pi/m)
            S += abs(c)**p; Sp += -p*abs(c)**(p-1)*(1 if c > 0 else -1)*s
        r = S**(-1/p); rp = -r*Sp/(p*S)
        return sqrt(r*r+rp*rp)
    return 2*m*quad(f, [0, pi/(2*m), pi/m])
for h in [mpf('0.004'), mpf('0.002'), mpf('0.001')]:
    d3 = (P(1+2*h)-2*P(1+h)+2*P(1-h)-P(1-2*h))/(2*h**3)
    print("h", h, " FD P'''(1) =", nstr(d3, 15))
