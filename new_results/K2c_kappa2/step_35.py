work='<project>/research/08_theory/K2c_kappa2/'
code = open(work+'step_24.py').read().split("angles =")[0]
code += r'''
a = mp.pi/2
A = mp.sqrt(2*a*mp.sin(a))
b1 = beta1_formula(a); b2f = beta2_formula(a)
print('N     sig      Re[(P-1-b1 s-b2 s^2)/s^3]   Im')
for N in [8000, 12000, 18000, 26000, 38000, 56000]:
    th = mp.sqrt(2*a*N)*(1+1j)
    P = -(2/A)*N*mp.e**th*Iplus(N, a)
    sig = 1/mp.sqrt(mpf(N))
    b3eff = (P-1-b1*sig-b2f*sig**2)/sig**3
    print('%6d %.5f  %s' % (N, float(sig), mp.nstr(b3eff, 18)))
'''
open(work+'step_26.py','w').write(code)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python', work+'step_26.py'],capture_output=True,text=True,timeout=3600)
print(r.stdout); print(r.stderr[-800:])
