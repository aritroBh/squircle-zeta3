code = '''"""RST1R step_02: independent verification of R1 (modulus-exponent lemma).
(a) Pade log(1+u) >= 2u/(2+u): sympy derivative identity.
(b) rho <= r+1/2.
(c) |s+q| >= rho + q/sqrt2: exact difference identity + arb stress.
(d) (5): exact exponent <= Etilde(r), arb stress (all logs via bound methods).
"""
from flint import arb, acb, ctx
import sympy as sp
ctx.dps = 60

u = sp.symbols('u', nonnegative=True)
g = sp.log(1+u) - 2*u/(2+u)
gp = sp.simplify(sp.diff(g, u))
print("Pade g'(u) =", gp)
assert sp.simplify(gp - u**2/((1+u)*(2+u)**2)) == 0
print("Pade: g(0)=0, g'(u)=u^2/((1+u)(2+u)^2)>=0 => log(1+u)>=2u/(2+u) for u>=0. PROVED")

PI = arb.pi(); SQ2 = arb(2).sqrt(); a = PI/(2*SQ2)
W = acb(1,1)/SQ2

rs = [arb(0), arb('0.01'), arb('0.5'), arb(1), arb(2), arb(5), arb(10), arb(30), arb(100), arb(1000), arb('1e6')]
qs = [arb(1), arb(2), arb(3), arb(13), arb(50), arb(100)]
js = [0, 1, 2, 50, 99, 400, 10**5]
bad = 0; worst_d = None; min_c_abs = None
for r in rs:
    s = arb(1)/2 + r*W
    rho = abs(s)
    assert rho <= r + arb(1)/2 + arb('1e-50')
    for q in qs:
        diff = abs(s+q)**2 - (rho+q/SQ2)**2
        alg = q + q**2/2 + SQ2*q*(r-rho)
        assert abs(diff-alg) < arb('1e-35')*(1+abs(diff)), (r,q,diff,alg)
        low = q*(1-1/SQ2)+q**2/2
        assert diff > low - arb('1e-35'), (r,q)
        min_c_abs = diff if min_c_abs is None else min(min_c_abs, diff)
        for j in js:
            if j > 0:
                exact = -a*r + j*(rho/(rho+q/SQ2)).log()
            else:
                exact = -a*r
            Etil  = -a*r - SQ2*j*q/(2*r+1+q/SQ2)
            marg = Etil - exact
            if not (marg >= 0):
                bad += 1; print("VIOLATION", r.str(4), q.str(4), j, marg.str(6))
            worst_d = marg if worst_d is None else min(worst_d, marg)
print("(b) rho <= r+1/2: OK on grid (11 radii)")
print(f"(c) |s+q| >= rho+q/sqrt2: identity diff = q+q^2/2+sqrt2 q (r-rho) verified exactly;")
print(f"    lower bound q(1-1/sqrt2)+q^2/2 OK; min diff on grid = {min_c_abs.str(8)} > 0")
print(f"(d) (5): min margin Etilde - exact on grid = {worst_d.str(8)}; violations = {bad}")
assert bad == 0
print("R1 numeric stress: PASS")
'''
open('step_02.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run([PY,'step_02.py'],capture_output=True,text=True)
print("rc:",r.returncode); print(r.stdout); print(r.stderr[-2000:])
