from mpmath import mp, mpf, pi, sin, cos, tan, log, quad, polylog
mp.dps = 60
import sys; sys.path.insert(0,'.')
from squircle import H1u, L1

def Cl2(theta):
    return (polylog(2, mp.e**(1j*theta))).imag
def L1_closed(phi):
    return mpf('0.5')*(tan(phi)*log(sin(phi)) + (Cl2(pi-2*phi)+Cl2(2*phi))/2) - phi/2

for phi in [0.3, 0.8, 1.2, 1.5]:
    a = L1(phi); b = L1_closed(phi)
    print(phi, mp.nstr(a-b, 5))
