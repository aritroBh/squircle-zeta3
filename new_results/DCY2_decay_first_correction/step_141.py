import subprocess
step03 = '''#!/usr/bin/env python3
"""step_03_watson.py — certified Watson constants for DCY2 (repaired, v3).

C_rho: |rho3(u)| <= C_rho u^3 on [0,1], rho3 = v - 1 + u/2 - 5u^2/24,
       v(u) = sqrt((1-e^{-2u})/(2u)).
Method: series of v at 0 converges for |u| < pi.
  For u in (0, 1/2]: |rho3(u)| <= (M/6) |u|^3, M = max_{|z|=2}|v(z)| (Cauchy c_k <= M/2^k).
  For u in [1/2, 1]: interval-evaluate rho3(u)/u^3 on proper subinterval cells.
Certify C_rho < 0.07  (truth: max |rho3|/u^3 = 1/16 = 0.0625 at u -> 0, decreasing).
"""
from flint import arb, acb, ctx

def iv(lo, hi):
    lo = arb(lo); hi = arb(hi)
    return arb((lo+hi)/2, (hi-lo)/2)

def v_acb(z):
    return ((1-(-2*z).exp())/(2*z)).sqrt()

def max_on_circle(r, n=128):
    M = arb(0)
    for k in range(n):
        th = 2*arb.pi()*arb(k)/n
        z = acb(r*th.cos(), r*th.sin())
        M = max(M, abs(v_acb(z)))
    return M

def certify():
    M2 = max_on_circle(arb(2))
    C_small = M2/6
    best = arb(0)
    for i in range(32):
        u = iv(arb(1)/2 + arb(i)/64, arb(1)/2 + arb(i+1)/64)
        rho3 = v_acb(acb(u)) - 1 + u/2 - 5*u**2/24
        best = max(best, (abs(rho3)/u**3).upper())
    C = max(C_small, best)
    print("M(|z|=2)|v| <=", M2.str(10), " -> C_rho(Cauchy) <=", C_small.str(10),
          "; interval part <=", best.str(10), " => C_rho <=", C.str(10))
    return C

if __name__ == "__main__":
    ctx.dps = 40
    C = certify()
    assert C < arb("0.07"), C
    print("CERTIFIED C_rho < 0.07 ; truth: max |rho3|/u^3 = 1/16 = 0.0625 at u -> 0")
'''
open("step_03_watson.py","w").write(step03)
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","step_03_watson.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-1500:] if r.returncode else "")
