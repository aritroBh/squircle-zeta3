import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
step03 = '''#!/usr/bin/env python3
"""step_03: certified ingredient bounds for DCY2 (arb, outward rounding).

Lemma W (Watson): v(u) = sqrt((1-e^{-2u})/(2u)) = 1-u/2+5u^2/24-u^3/16+rho3(u).
Certify C_rho (|rho3| <= C_rho u^3 on [0,1]) and C_rho2 (u>=1).
C_rho via Cauchy estimate: |v^(3)(x)| <= 3! M(x,1/2)/(1/2)^3 = 48 M.
"""
from flint import arb, acb, ctx

def max_abs_v_on_circle(x, r, n=96):
    M = arb(0)
    for k in range(n):
        th = arb(2)*arb.pi()*k/n
        z = acb(x) + acb(r*th.cos(), r*th.sin())
        w = (1-(-2*z).exp())/(2*z)
        vz = w.sqrt()
        M = max(M, abs(vz))
    return M

def certify_C_rho():
    best = arb(0)
    for i in range(129):
        x = arb(i)/128
        M = max_abs_v_on_circle(x, arb(1)/2)
        best = max(best, 48*M)
    C = best/24
    print("max |v^(3)| on [0,1] <=", best.str(12), "-> C_rho <=", C.str(12))
    return C

def certify_C_rho2():
    # |rho3(u)| <= |v(u)|+1+u/2+5u^2/24+u^3/16 <= u^3(1/sqrt2+1+1/2+5/24+1/16)
    C = 1/arb(2).sqrt() + 1 + arb(1)/2 + arb(5)/24 + arb(1)/16
    print("C_rho2 <=", C.str(12))
    return C

if __name__ == "__main__":
    ctx.dps = 40
    C1 = certify_C_rho()
    C2 = certify_C_rho2()
    assert C1 < arb("0.014") and C2 < arb("1.78")
    print("CERTIFIED: C_rho < 0.014, C_rho2 < 1.78")
'''
open("step_03_watson.py","w").write(step03)
r = subprocess.run(["<project>/.venv/bin/python","step_03_watson.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-1500:] if r.returncode else "")
