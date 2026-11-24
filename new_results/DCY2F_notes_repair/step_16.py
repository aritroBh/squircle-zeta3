import sys, time
sys.path.insert(0, "<project>/research/08_theory/DCY2_decay_first_correction")
from flint import arb, acb, ctx
ctx.dps = 35
t0=time.monotonic()
from step_06_finite import main_arb, residual_check, varpi1
PI = arb.pi()
N = arb(381); th = (PI*N).sqrt()
M = main_arb(N)
print("M (arb)      =", M.str(15))
# their T1
T1 = -(PI.sqrt()/N)*(-th).exp()*( (th).cos()*(1+varpi1.real/N.sqrt()) - (th).sin()*(varpi1.imag/N.sqrt()) )
print("their T1     =", T1.str(15))
X = (PI.sqrt()/N)*(-th).exp()*( (th).cos()*(1+varpi1.real/N.sqrt()) - (th).sin()*(varpi1.imag/N.sqrt()) )
print("X = -T1      =", X.str(15))
scale = N**2*th.exp()
print("res=|M-T1| scaled =", (abs(M-T1)*scale).str(10))
print("|M+X| scaled      =", (abs(M+X)*scale).str(10))
print("|M-X| scaled      =", (abs(M-X)*scale).str(10))
print("elapsed", time.monotonic()-t0)
