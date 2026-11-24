import subprocess, os, time
wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
code = r'''
import flint
from flint import acb, arb
import mpmath as mp
mp.mp.dps = 50

def I_cert(eta_str, U=64, tol='1e-35'):
    eta = acb(eta_str)
    ctx = flint.ctx
    ctx.prec = 160
    def f(u, ana):
        return (1+u.exp())**(-2-eta) * 2 * (-eta*u/2).exp() * (eta*u).cosh().sqrt() * (-u).exp()
    ctx.dps = 50
    v, err = acb.calc_integrate(f, acb(0), acb(U), rel_tol=acb(tol), abs_tol=acb(tol))
    # rigorous tail bound on [U, inf): |integrand| <= 2 * 2^{2+|eta|} e^{-u}  (|H|<=2*2^{|eta|/... })
    # |H(eta u)| = |2 e^{-eta u/2} sqrt(cosh eta u)| <= 2 e^{|eta| u/2} e^{|eta|u/2} = 2 e^{|eta| u}
    # |K| = (1+e^{-u})^{|eta+2|-sign} <= max(1, 2^{|eta+2|})
    # tail <= 2*max(1,2^{|eta+2|}) * e^{-U} / (1 - |eta|)   (if |eta|<1)
    import math
    em = abs(float(eta_str))
    tail = 2*max(1.0, 2.0**(abs(2+em)))*math.exp(-U)/(1-em)
    return v, err, tail

for eta in ['0.1','0.3','-0.1','-0.3']:
    v, err, tail = I_cert(eta)
    print(f"eta={eta}: calc_integrate:")
    print("   value:", v)
    print("   abs err (calc):", err, " tail bound:", tail)
'''
with open(wd+"/step_30.py","w") as fp: fp.write(code)
t0=time.time()
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_30.py"], capture_output=True, text=True, cwd=wd, timeout=1700)
print(r.stdout[-3000:]); print("STDERR:", r.stderr[-1500:]); print("elapsed", time.time()-t0)
