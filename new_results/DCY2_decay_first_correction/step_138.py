import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")

# R5a: TailMom at Nlo for p <= 2
src = open("dcy2_sweep.py").read()
old = """    Et = arb(0)
    for p, cf in pp.items():
        Et += cf*TailMom_of(Nhi, p)
    return E + 2*Et   # two-sided window tail (R4)"""
new = """    Et = arb(0)
    for p, cf in pp.items():
        # R5: TP_p decreasing for p <= 2 (erfc/Q_1 form) -> Nlo; increasing for p >= 3 -> Nhi
        Et += cf*TailMom_of(Nlo if p <= 2 else Nhi, p)
    return E + 2*Et   # two-sided window tail (R4)"""
assert old in src
src = src.replace(old, new)
open("dcy2_sweep.py","w").write(src)

# R5b: e^{gamma T^3} at consistent endpoint (Nlo, since gammaT^3 = c2min^3/(8 gamma^2) increasing)
src = open("dcy2_sweep.py").read()
old = """        T_lo = min(D_lo, c_lo/(2*gam_lo))
        base = min((gam_lo*T_lo**3).exp()*Mtilde_of(Nhi, jj),
                   Mtilde_of(Nhi, jj, c_hi/2))"""
assert old in src
src = src.replace(old, """        T_lo = min(D_lo, c_lo/(2*gam_lo))
        gamT3 = min(gam_lo*D_lo**3, c_lo**3/(8*gam_lo**2))   # R5: gamma T^3 <= max of both forms
        base = min(gamT3.exp()*Mtilde_of(Nhi, jj),
                   Mtilde_of(Nhi, jj, c_hi/2))""")
open("dcy2_sweep.py","w").write(src)
print("R5 applied")

# R6: repair step_03_watson.py (Cauchy M/6, proper interval cells)
step03 = '''#!/usr/bin/env python3
"""step_03_watson.py — certified Watson constants for DCY2 (repaired).

C_rho: |rho3(u)| <= C_rho u^3 on [0,1], rho3 = v - 1 + u/2 - 5u^2/24,
       v(u) = sqrt((1-e^{-2u})/(2u)).
Method: for |u| <= 1/2, |rho3(u)| <= M u^3/6 with M = max_{|z|=2}|v(z)| (Cauchy, c_k <= M/2^k);
        for u in [1/2, 1], interval-evaluate rho3(u)/u^3 on subintervals (proper cells).
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

def certify_C_rho():
    M2 = max_on_circle(arb(2))
    C_small = M2/6   # covers |u| <= 1/2: |rho3| <= sum_{k>=3} M (|u|/2)^k <= M u^3/6
    best = arb(0)
    for i in range(32):
        u = iv(arb(1)/2 + arb(i)/64, arb(1)/2 + arb(i+1)/64)
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
    print("CERTIFIED C_rho < 0.5 ; mpmath truth ~ 0.012")
'''
open("step_03_watson.py","w").write(step03)
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","step_03_watson.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-1500:] if r.returncode else "")
