import sys, time
sys.path.insert(0, '../DCY2_decay_first_correction')
from mpmath import mp, mpf, pi, sqrt, exp, log
import dcy2_audit as A

def J_direct_hp(N, dps=60, seg_factor=40):
    mp.dps = dps
    from dcy2_common import Hc, omega
    j=1j
    rho = sqrt(2*N/pi)
    f = lambda x: Hc(omega*x)*exp(j*pi*omega*x/2)*(omega*x/(1+omega*x))**N/(omega*x)*omega
    pts = [0, rho/8, rho/2, rho, 2*rho, 4*rho, 10*rho, seg_factor*rho]
    v = mp.quad(f, pts)
    # tail beyond seg_factor*rho: |integrand| <= |H| e^{-bR x} -> tiny
    return (2/pi)*v

for N in [400, 2000, 10000]:
    mp.dps = 60
    tq = A.true_quant(N)
    model = tq['Model']
    for dp in [40, 60, 80]:
        J = J_direct_hp(N, dp)
        d = abs(J-model)*N**2*exp(sqrt(pi*N))
        print(f"N={N}, dps={dp}: |J-model|N^2e^th (direct quad) = {mp.nstr(d,8)}")
    print()
