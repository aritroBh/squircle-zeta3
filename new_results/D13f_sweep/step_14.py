
import subprocess, os
os.chdir("research/08_theory/D13f_sweep")
PY = "<home>)/.venv/bin/python"
test = r'''
from flint import arb, acb, ctx
ctx.dps = 70
pi = arb.pi(); rt2 = arb(2).sqrt(); a = pi/(2*rt2)
N = arb(100)
x0 = arb("7.92549")
x2 = x0/2
RePhi = -a*x2 + N*x2.log() - N*(1+rt2*x2+x2*x2).log()/2
th = (pi*N).sqrt()
vf = rt2/x2 + 1/x2**2
ReX2 = -a*x2 + th - (N/2)*(1+vf).log()
print("N log x2 =", N*x2.log(), "  th =", th)
print("step_34 RePhi(x2) =", RePhi)
print("step_36 ReX2      =", ReX2)
print("ReX2 - th =", ReX2 - th)
# so is N log x2 <= th ???
print("N log x2 <= th ?", arb(N*x2.log()) <= th)
# Hmm: maybe they intended RePhi(x2) = -a x2 + N log x2 - (N/2) log(...) 
# and bound N log x2 by th? That's wrong. Check what makes TLf ~4.7/sc:
lam2 = -a + N/x2 - N*(x2+1/rt2)/(1+rt2*x2+x2*x2)
sqpi = pi.sqrt()
def Hr(t): return sqpi/4*(t/2).gamma()/(t/2+arb("1.5")).gamma()
G1h = Hr(1/rt2)
sc = N**arb("-1.5")*(-th).exp()
TLf34 = G1h*RePhi.exp()/lam2
TLf36 = G1h*ReX2.exp()/lam2
print("TLf34/sc =", TLf34/sc)
print("TLf36/sc =", TLf36/sc)
'''
r = subprocess.run([PY, "-c", test], capture_output=True, text=True)
print(r.stdout)
print("STDERR:", r.stderr[:2000])
