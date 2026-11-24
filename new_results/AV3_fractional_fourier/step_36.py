from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, ci, zeta, euler, psi, pslq, quad, tan)
mp.dps = 100
def N(x,d=42): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
Lp = log(pi); eta3 = 3*zeta(3)/4; zp3 = zeta(3, derivative=1)
def bp2(n):
    Jc = -( (ci(n*pi)-euler-((-1)**n)*Lp-log(n))/n + (1-(-1)**n)/n )/n
    Kc = ((-1)**n - 1)/n**2
    return (2/pi)*( -pi**2*Lp*(-1)**n/n + (2*Jc + Kc)/n )
CiS = nsum(lambda n: ci(n*pi)/n**3, [1, inf])
Wp2_direct = nsum(lambda n: bp2(n), [1, inf])
Wp2_pieces = 2*pi*Lp*log(2) - (4/pi)*(CiS + zp3 + Lp*eta3 - euler*zeta(3)) - (21/(2*pi))*zeta(3)
Wp2_int = (1/pi)*quad(lambda t: t**2*log(t)/tan(t/2), [0, pi])
print("W'(2) direct nsum =", N(Wp2_direct,45), flush=True)
print("W'(2) pieces      =", N(Wp2_pieces,45), flush=True)
print("W'(2) integral    =", N(Wp2_int,45), flush=True)
basis = [mpf(1), Lp*log(2)*pi, pi*log(2), pi*log(2)**2, zeta(3)/pi, zeta(3)*Lp/pi,
         euler*zeta(3)/pi, zeta(3)/pi*log(2), zp3/pi, zeta(5)/pi, pi**3*log(2), pi]
names = ['1','pi Lp log2','pi log2','pi log2^2','z3/pi','z3 Lp/pi','gam z3/pi','z3 log2/pi',"z3'/pi",'z5/pi','pi^3 log2','pi']
rel = pslq([Wp2_direct]+basis, tol=mpf('1e-95'), maxcoeff=10**12, maxsteps=10000)
print("PSLQ W'(2):", rel, flush=True)
if rel:
    for c,nm in zip(rel, ['target']+names):
        if c: print("  ", c, nm, flush=True)
basis2 = [mpf(1), log(2), Lp, euler, pi**2, pi**3, zeta(3), zeta(3)*log(2), zeta(3)*Lp, euler*zeta(3), zp3, zeta(3)**2, pi**3*log(2), pi**3*Lp]
names2 = ['1','log2','Lp','gam','pi^2','pi^3','z3','z3 log2','z3 Lp','gam z3',"z3'",'z3^2','pi^3 log2','pi^3 Lp']
rel2 = pslq([CiS]+basis2, tol=mpf('1e-95'), maxcoeff=10**12, maxsteps=10000)
print("PSLQ CiS:", rel2, flush=True)
if rel2:
    for c,nm in zip(rel2, ['target']+names2):
        if c: print("  ", c, nm, flush=True)
logline("step_32: W'(2) three methods; PSLQ on W'(2) and CiS")
