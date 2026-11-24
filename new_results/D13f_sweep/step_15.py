
import subprocess, os
os.chdir("research/08_theory/D13f_sweep")
PY = "<home>)/.venv/bin/python"
test = r'''
from mpmath import mp, mpf, mpc, pi, sqrt, exp, log, gamma
mp.dps = 50
b = exp(mpc(0,pi/4)); sqpi = sqrt(pi)
def H(s): return sqpi/4*gamma(s/2)/gamma(s/2+mpf("1.5"))
def Hr(t): return sqpi/4*gamma(t/2)/gamma(t/2+mpf("1.5"))
def g(x): return H(b*x)/x
for x in [7.9, 5.3, 10.0, 3.0, 1.0, 40.0]:
    lhs = abs(g(x))
    rhs = Hr(x/sqrt(2))/x
    print(f"x={x}: |g(x)|={mp.nstr(lhs,12)}  Hr(x/rt2)/x={mp.nstr(rhs,12)}  rhs>=lhs: {rhs>=lhs}  ratio={mp.nstr(lhs/rhs,6)}")
'''
r = subprocess.run([PY, "-c", test], capture_output=True, text=True)
print(r.stdout)
print("STDERR:", r.stderr[:2000])
