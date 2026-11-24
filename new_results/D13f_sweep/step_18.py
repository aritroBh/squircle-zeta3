
import subprocess, os
os.chdir("research/08_theory/D13f_sweep")
PY = "<home>)/.venv/bin/python"
test = r'''
from mpmath import mp, mpf, mpc, pi, sqrt, exp, log, gamma
mp.dps = 50
b = exp(mpc(0,pi/4)); rt2 = sqrt(2)
def gap(x):   # gamma'(x)
    s=b*x
    return b*(mp.polygamma(0,s/2)-mp.polygamma(0,s/2+mpf("1.5")))/2 - 1/x
def gapp(x):  # gamma''(x)
    s=b*x
    return b*b*(mp.polygamma(1,s/2)-mp.polygamma(1,s/2+mpf("1.5")))/4 + 1/x**2
def psi1(t): return mp.polygamma(1,t)   # trigamma
def psi2(t): return mp.polygamma(2,t)   # tetragamma (negative)
print("x      |gap|    (3/4)psi'(t)+1/x   (3/4)(1/t^2+1/t)+1/x    |gapp|   (3/8)|psi''|+1/x^2   (3/8)(2/t^3+1/t^2)+1/x^2")
for x in [0.5, 1.0, 2.0, 4.0, 7.9, 16.0, 40.0, 200.0]:
    t = x/(2*rt2)
    g1 = abs(gap(x)); b1 = 0.75*psi1(t)+1/x; b1e = 0.75*(1/t**2+1/t)+1/x
    g2 = abs(gapp(x)); b2 = 0.375*abs(psi2(t))+1/x**2; b2e = 0.375*(2/t**3+1/t**2)+1/x**2
    print(f"{x:6} {mp.nstr(g1,6):>10} {mp.nstr(b1,6):>15} {mp.nstr(b1e,6):>15}   {mp.nstr(g2,6):>10} {mp.nstr(b2,6):>15} {mp.nstr(b2e,6):>15}")
# Gh monotonicity
def Gh(eta, rho, y):
    t = rho*y
    return rho*(eta+t*rt2)/(rt2*eta*(eta*eta+rt2*eta*t+t*t))
print("Gh values (rho~1, y=0.125):", [mp.nstr(Gh(e,1.0,0.125),6) for e in [0.4,0.6,0.8,1.0,1.2,1.5,2.0,3.0]])
'''
r = subprocess.run([PY, "-c", test], capture_output=True, text=True)
print(r.stdout)
print("STDERR:", r.stderr[:1500])
