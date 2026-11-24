from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, ci, zeta, euler, psi, pslq)
mp.dps = 120
def N(x,d=42): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
Lp = log(pi); eta3 = 3*zeta(3)/4; zp3 = zeta(3, derivative=1)
CiS = nsum(lambda n: ci(n*pi)/n**3, [1, inf])
Wp2 = 2*pi*Lp*log(2) - (4/pi)*(CiS + zp3 + Lp*eta3 - euler*zeta(3)) - (21/(2*pi))*zeta(3)
print("W'(2) =", N(Wp2,50), flush=True)
print("CiS   =", N(CiS,50), flush=True)

def try_pslq(target, basis, names, tag):
    x = [target]+basis
    rel = pslq(x, tol=mpf('1e-100'), maxcoeff=10**10, maxsteps=10000)
    print(tag, "rel:", rel, flush=True)
    if rel:
        resid = abs(sum(c*v for c,v in zip(rel,x)))
        norm  = max(abs(c) for c in rel)
        print("  residual:", nstr(resid,3), " maxcoeff:", nstr(norm,3), flush=True)
        for c,nm in zip(rel, ['target']+names):
            if c: print("   ", c, nm, flush=True)

basis = [mpf(1), Lp*log(2)*pi, pi*log(2), pi*log(2)**2, zeta(3)/pi, zeta(3)*Lp/pi,
         euler*zeta(3)/pi, zeta(3)/pi*log(2), zp3/pi, zeta(5)/pi, pi**3*log(2), pi, zeta(3)**2/pi]
names = ['1','pi Lp log2','pi log2','pi log2^2','z3/pi','z3 Lp/pi','gam z3/pi','z3 log2/pi',"z3'/pi",'z5/pi','pi^3 log2','pi','z3^2/pi']
try_pslq(Wp2, basis, names, "W'(2)")
basis2 = [mpf(1), log(2), Lp, euler, pi**2, pi**3, zeta(3), zeta(3)*log(2), zeta(3)*Lp, euler*zeta(3), zp3, zeta(3)**2, pi**3*log(2), pi**3*Lp]
names2 = ['1','log2','Lp','gam','pi^2','pi^3','z3','z3 log2','z3 Lp','gam z3',"z3'",'z3^2','pi^3 log2','pi^3 Lp']
try_pslq(CiS, basis2, names2, "CiS")
logline("step_33: accurate W'(2), CiS; PSLQ with residual norms")
