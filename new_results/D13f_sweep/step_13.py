
import subprocess, os
os.chdir("research/08_theory/D13f_sweep")
PY = "<home>)/.venv/bin/python"
test = r'''
from flint import arb, acb, ctx
ctx.dps = 70
pi = arb.pi(); rt2 = arb(2).sqrt()
b = (acb(0,1)*pi/4).exp()
def iv(lo, hi):
    lo = arb(lo).lower(); hi = arb(hi).upper()
    return arb((lo+hi)/2, (hi-lo)/2)
# x0 box like dbg1: [7.9 +/- 0.0789]
x0 = iv(7.82, 7.98)
print("x0:", x0)
s = b*x0
print("s:", s)
z1 = s/2
z2 = s/2 + arb("1.5")
print("z1:", z1)
d1 = z1.digamma()
d2 = z2.digamma()
print("digamma(z1):", d1)
print("digamma(z2):", d2)
diff = d1 - d2
print("diff:", diff)
gap = b*diff/2 - 1/x0
print("gap:", gap, " abs:", abs(gap))
# compare: individual widths
print("width digamma(z1):", d1.real.rad(), d1.imag.rad())
# lgamma similarly
lg1 = z1.lgamma(); lg2 = z2.lgamma()
print("lgamma(z1):", lg1)
print("lgamma diff:", lg1-lg2)
'''
r = subprocess.run([PY, "-c", test], capture_output=True, text=True)
print(r.stdout)
print("STDERR:", r.stderr[:2000])
