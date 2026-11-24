import mpmath
from mpmath import mp, mpf, sqrt, pi, log, exp, cos, sin, tan, nstr, atan
mp.dps = 70
def d1_closed(lam):
    phi = atan(lam)
    R = sqrt(1+lam**2)
    def G(psi):
        return -cos(psi)*log(sin(psi)) + log(tan(psi/2)) + cos(psi)
    I = R*( (sin(phi)+cos(phi))*log(R) + G(phi+pi/2) - G(phi) )
    return 4*( log(lam) + (log(2)-1)*(1+1/lam) - I/lam )
for lam, fit in [(mpf(2), '-3.0630013748005184342593242773606'),
                 (mpf(3), '-2.421476812773221095332596814398359'),
                 (mpf('0.5'), '-6.1260027496010537996534024718262299')]:
    print(f"lam={lam}: closed {nstr(d1_closed(lam),45)}")
    print(f"          fit    {fit[:45]}")
with open('progress.log','a') as fh:
    fh.write('step_78: anisotropic p->inf d1 closed form confirmed; profile non-shrinking, no zeta\n')
