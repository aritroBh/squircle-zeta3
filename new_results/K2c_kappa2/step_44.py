work='<project>/research/08_theory/K2c_kappa2/'
code = open(work+'step_30.py').read().split("# model pieces")[0]
code += r'''
import numpy as np
# d1 from symbolic formula (step_32): d1 = (q^2*(-9 sqrt2 g1 + 18 g1 - q^2 + 2 sqrt2 q^2 - 15 sqrt2 + 15) - 3 sqrt2 + 6)/(12 q), tot = -(1+d1 sig)
def d1_formula(Z, g1):
    q = mp.sqrt(Z); sq2 = mp.sqrt(2)
    return -( (q**2*(-9*sq2*g1 + 18*g1 - q**2 + 2*sq2*q**2 - 15*sq2 + 15) - 3*sq2 + 6)/(12*q) )
    # minus because totn = -(1+d1 sig): C_img = -e^{Z/2}(1 + d1 sig) with d1 = -d1_printed

for a in [mp.pi/3, mp.pi/2, 2*mp.pi/3, 3*mp.pi/4, 5*mp.pi/6]:
    Z = 1j*a; g1 = g1f(a)
    rows = []
    for N in [1500, 3000, 6000, 12000]:
        sig = 1/mp.sqrt(mpf(N))
        th = mp.sqrt(2*a*N)
        I1 = Iplus(N, a, 1); I2 = Iplus(N, a, 2)
        C = -2*I2/I1*exp((s2-1)*th*(1+1j))
        rows.append((sig, C/(-exp(Z/2))))
    Arows = [[s, s*s] for s,_ in rows]
    yv = [complex(v-1) for _,v in rows]
    sol, *_ = np.linalg.lstsq(np.array(Arows,dtype=complex), np.array(yv,dtype=complex), rcond=None)
    d1ext = mpc(sol[0].real, sol[0].imag)
    d1f = d1_formula(Z, g1)
    print('alpha=%s' % mp.nstr(a,5))
    print('  d1 extracted = %s   formula = %s   |diff|=%s' % (mp.nstr(d1ext,16), mp.nstr(d1f,16), mp.nstr(abs(d1ext-d1f),3)))
'''
open(work+'step_33.py','w').write(code)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python', work+'step_33.py'],capture_output=True,text=True,timeout=5400)
print(r.stdout); print(r.stderr[-500:])
