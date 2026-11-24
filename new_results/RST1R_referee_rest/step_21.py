code = '''"""RST1R step_02: independent verification of R1 (modulus-exponent lemma).
(a) Pade log(1+u) >= 2u/(2+u): sympy derivative + arb stress.
(b) rho^2 = r^2 + r/sqrt2 + 1/4; rho <= r+1/2.
(c) |s+q| >= rho + q/sqrt2: exact difference identity + arb stress (worst case large r).
(d) (5): -a r + j log(rho/(rho+q/sqrt2)) <= -a r - sqrt2 j q/(2r+1+q/sqrt2), arb stress.
"""
from flint import arb, acb, ctx
import sympy as sp
ctx.dps = 60

# (a) Pade via sympy: g(u)=log(1+u)-2u/(2+u), g(0)=0, g'(u)=u^2/((1+u)(2+u)^2) >= 0
u = sp.symbols('u', nonnegative=True)
g = sp.log(1+u) - 2*u/(2+u)
gp = sp.simplify(sp.diff(g, u))
print("Pade g'(u) =", gp)
assert sp.simplify(gp - u**2/((1+u)*(2+u)**2)) == 0
print("Pade: g(0)=0, g' = u^2/((1+u)(2+u)^2) >= 0  =>  log(1+u) >= 2u/(2+u). PROVED")

PI = arb.pi(); SQ2 = arb(2).sqrt(); a = PI/(2*SQ2)
W = acb(1,1)/SQ2

# (b)(c)(d) arb stress
worst_c = arb(10**9)   # min margin of |s+q|^2 - (rho+q/sqrt2)^2 relative
worst_d = arb(10**9)   # min of Etilde - exact
min_c_abs = arb(10**9)
import itertools
rs = [arb(0), arb('0.01'), arb('0.5'), arb(1), arb(2), arb(5), arb(10), arb(30), arb(100), arb(1000), arb('1e6')]
qs = [arb(1), arb(2), arb(3), arb(13), arb(50), arb(100)]
js = [0, 1, 2, 50, 99, 400, 10**5]
bad = 0
for r in rs:
    s = arb(1)/2 + r*W
    rho = abs(s)
    # (b) rho <= r+1/2
    assert rho <= r + arb(1)/2 + arb('1e-50')
    for q in qs:
        lhs2 = abs(s+q)**2; rhs2 = (rho+q/SQ2)**2
        diff = lhs2 - rhs2
        # algebra: diff = q + q^2/2 + sqrt2 q (r - rho)
        alg = q + q**2/2 + SQ2*q*(r-rho)
        assert abs(diff-alg) < arb('1e-40')*(1+abs(diff)), (r,q,diff,alg)
        low = q*(1-1/SQ2)+q**2/2   # claimed lower bound using r-rho >= -1/2
        assert diff > low - arb('1e-40'), (r,q)
        if diff > 0: min_c_abs = min(min_c_abs, diff)
        for j in js:
            exact = -a*r + j*acb.log(rho/(rho+q/SQ2)).real if j>0 else -a*r
            Etil  = -a*r - SQ2*j*q/(2*r+1+q/SQ2)
            marg = Etil - exact
            if marg < 0: bad += 1; print("VIOLATION", r,q,j,marg)
            worst_d = min(worst_d, marg)
print(f"(b) rho <= r+1/2: OK on grid")
print(f"(c) |s+q| >= rho+q/sqrt2: exact identity diff=q+q^2/2+sqrt2 q(r-rho) verified; min diff on grid = {min_c_abs.str(8)} > 0")
print(f"(d) (5) exact <= Etilde: min margin on grid = {worst_d.str(8)}; violations = {bad}")
assert bad == 0
print("R1 numeric stress: PASS")
'''
open('step_02.py','w').write(code)
import subprocess
r=subprocess.run(['<project>/.venv/bin/python','step_02.py'],capture_output=True,text=True)
print(r.stdout); print(r.stderr[-2000:] if r.returncode else "exit 0")
