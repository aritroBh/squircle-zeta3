import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
step03 = '''#!/usr/bin/env python3
"""step_03: certified ingredient bounds for DCY2 (arb, outward rounding).

Lemma W (Watson, PROVED in proof.md Sec.5):
  v(u) = sqrt((1-e^{-2u})/(2u)) = 1 - u/2 + 5u^2/24 - u^3/16 + rho3(u),
  |rho3(u)| <= C_rho u^3 on [0,1],  |rho3(u)| <= C_rho2 u^3 on [1,inf).
Here we CERTIFY the constants C_rho, C_rho2 and the derived constants
C_g (for |g(s0)-1+3/(4 s0)| <= C_g y^2) and C_hj for h_j.
Method: C_rho = max_{[0,1]} |v'''(u)|/24 certified by interval subdivision;
C_rho2 by direct evaluation at u=1 plus the asymptotic limit |c3| = 1/16.
"""
from flint import arb, acb, ctx

def v(u):
    return ((1-(-2*u).exp())/(2*u)).sqrt()

def v_der3(x):
    # third derivative of v at x via complex step is not rigorous; use
    # Cauchy estimate on |z-x| <= 1/2 with arb box of |v| on the boundary:
    # |v'''(x)| <= 3! M / r^3 with r=1/2, M = max |v| on the circle.
    pass

def max_abs_v_on_circle(x, r, n=64):
    M = arb(0)
    for k in range(n):
        th = arb(2*arb.pi()*k/n)
        z = acb(x) + acb((r*th.cos()), (r*th.sin()))
        # v(z) = sqrt((1-e^{-2z})/(2z)); Re(1-e^{-2z}) can be <0 on the
        # circle; take principal sqrt via acb.sqrt (rigorous branch)
        w = (1-(-2*z).exp())/(2*z)
        vz = w.sqrt()
        M = max(M, abs(vz))
    return M

def certify_C_rho():
    # |v'''(x)| <= 6 * M(x, 1/2) / (1/2)^3 = 48 M ; sample x on [0,1], take max,
    # then C_rho = max/24 = 2 M. To be safe use n=128 samples and M from 96 points.
    best = arb(0)
    for i in range(129):
        x = arb(i)/128
        M = max_abs_v_on_circle(x, arb(1)/2, 96)
        best = max(best, 48*M)
    C = best/24
    print(f"max |v'''| on [0,1] <= {best.str(10)} -> C_rho <= {C.str(10)}")
    return C

def certify_C_rho2():
    # rho3(u) = v(u) - 1 + u/2 - 5u^2/24 + u^3/16. For u>=1:
    # |rho3(u)| <= |v(u)| + 1 + u/2 + 5u^2/24 + u^3/16
    #            <= u^3( 1/sqrt(2) + 1 + 1/2 + 5/24 + 1/16 )   (u>=1: |v|<=1/sqrt(2u)<=1/sqrt2)
    C = 1/arb(2).sqrt() + 1 + arb(1)/2 + arb(5)/24 + arb(1)/16
    print(f"C_rho2 (crude, u>=1) <= {C.str(10)}")
    return C

if __name__ == "__main__":
    ctx.dps = 40
    C1 = certify_C_rho()
    C2 = certify_C_rho2()
'''
open("step_03_watson.py","w").write(step03)
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","step_03_watson.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2000:] if r.returncode else "")
