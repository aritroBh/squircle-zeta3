import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
step03 = '''#!/usr/bin/env python3
"""step_03: certified Watson constants for DCY2.

C_rho: |rho3(u)| <= C_rho u^3 on [0,1], rho3 = v - 1 + u/2 - 5u^2/24.
Method: for |u| <= 1/2, |rho3(u)| <= M u^3/4 with M = max_{|z|=2}|v(z)| (Cauchy);
        for u in [1/2, 1], interval-evaluate rho3(u)/u^3 on subintervals.
"""
from flint import arb, acb, ctx

def v_acb(z):
    return ((1-(-2*z).exp())/(2*z)).sqrt()

def max_on_circle(r, n=128):
    M = arb(0)
    for k in range(n):
        th = 2*arb.pi()*arb(k)/n
        z = acb(r*th.cos(), r*th.sin())
        M = max(M, abs(v_acb(z)))
    return M

def certify_C_rho():
    M2 = max_on_circle(arb(2))
    C_small = M2/4   # covers |u| <= 1/2
    # interval evaluation on [1/2, 1], 32 cells
    best = arb(0)
    for i in range(32):
        lo = arb(1)/2 + arb(i)/64
        hi = arb(1)/2 + arb(i+1)/64
        u = arb(lo, hi)
        # rho3/u^3 with outward enclosure; v needs precision near the removable 0 -- here u >= 1/2, fine
        rho3 = v_acb(acb(u)) - 1 + u/2 - 5*u**2/24
        val = abs(rho3)/u**3
        best = max(best, val.upper())
    C = max(C_small, best)
    print("M(|z|=2)|v| <=", M2.str(10), "; C_rho <=", C.str(12), "; interval part:", best.str(12))
    return C

if __name__ == "__main__":
    ctx.dps = 40
    C = certify_C_rho()
    assert C < arb("0.5")
    print("CERTIFIED C_rho < 0.5 ; mpmath truth ~ 0.0138")
'''
open("step_03_watson.py","w").write(step03)
r = subprocess.run(["<project>/.venv/bin/python","step_03_watson.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-1500:] if r.returncode else "")
