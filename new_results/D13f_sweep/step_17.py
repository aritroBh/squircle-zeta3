
import subprocess, os
os.chdir("research/08_theory/D13f_sweep")
PY = "<home>)/.venv/bin/python"
test = r'''
from mpmath import mp, mpf, mpc, pi, sqrt, exp, log, gamma
mp.dps = 50
b = exp(mpc(0,pi/4)); a = pi/(2*sqrt(2)); rt2 = sqrt(2); sqpi = sqrt(pi)
def ReP1(x,N):  return -a + N/x - N*(x+1/rt2)/(1+rt2*x+x*x)
def ReP2(x,N):  return -N/x**2 + N*x*(x+rt2)/(1+rt2*x+x*x)**2
def x0_of(N): return mp.findroot(lambda x: ReP1(x,N), sqrt(2*N/pi))
def Gh(eta, rho, y):
    t = rho*y
    return rho*(eta+t*rt2)/(rt2*eta*(eta*eta+rt2*eta*t+t*t))
def lamR_reg(eta, rho, y):   # = -ReP1(eta*x0, N)
    return a - (pi*rho/2)*Gh(eta, rho, y)
def curv2(y, rho, eta):
    num = rt2*eta**3*rho + 4*eta**2*rho**2*y + 2*rt2*eta*rho**3*y**2 + rho**4*y**3
    den = eta**2*(eta**2+rt2*eta*rho*y+rho**2*y**2)**2
    return (pi*rho**2*y/4)*num/den
print("=== verify lamR = a - (pi rho/2) Gh(eta)  vs  -ReP1(eta x0,N) ===")
for Nv in [100, 400, 10000, 10**8]:
    N = mpf(Nv); x0 = x0_of(N); y0 = 1/x0; y = sqrt(pi/(2*N)); rho = y0/y
    for eta in [0.5, 0.9, 1.0, 1.1, 1.3, 2.0]:
        lhs = -ReP1(eta*x0, N)
        rhs = lamR_reg(eta, rho, y)
        print(f"  N={Nv} eta={eta}: -ReP1={mp.nstr(lhs,10)}  reg={mp.nstr(rhs,10)}  diff={mp.nstr(lhs-rhs,3)}")
print("=== verify curv2 vs |ReP2|/2 ===")
for Nv in [100, 10000]:
    N = mpf(Nv); x0 = x0_of(N); y0=1/x0; y = sqrt(pi/(2*N)); rho = y0/y
    for eta in [0.7, 1.0, 1.3]:
        lhs = abs(ReP2(eta*x0,N))/2
        rhs = curv2(y, rho, eta)
        print(f"  N={Nv} eta={eta}: |ReP2|/2={mp.nstr(lhs,10)}  curv2={mp.nstr(rhs,10)}  ratio={mp.nstr(lhs/rhs,6)}")
'''
r = subprocess.run([PY, "-c", test], capture_output=True, text=True)
print(r.stdout)
print("STDERR:", r.stderr[:1500])
