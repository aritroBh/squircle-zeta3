from mpmath import mp, mpf, zeta, log, pi, nprint, polylog, clsin, catalan, pslq, quad, sqrt
mp.dps = 100
P1 = quad(lambda t: log(4/(1-t**2))**2/(3+t**2), [0,1])
d1 = mpf(open('d1.txt').read())
# log-sine candidate constants
Ls3_60 = quad(lambda x: log(2*mp.sin(x/2))**2, [0, pi/3])
Ls3_90 = quad(lambda x: log(2*mp.sin(x/2))**2, [0, pi/2])
Ls3_120 = quad(lambda x: log(2*mp.sin(x/2))**2, [0, 2*pi/3])
print("Ls3(pi/3) =", nprint(Ls3_60, 30), " 7pi^3/108 =", nprint(7*pi**3/108, 30))
print("Ls3(pi/2) =", nprint(Ls3_90, 30))
print("Ls3(2pi/3)=", nprint(Ls3_120, 30))
G = catalan; z3 = zeta(3); L2 = log(2); C13 = clsin(2,pi/3); C16 = clsin(2,pi/6)
Li3h = polylog(3, mpf('0.5'))
basis = [1, pi, pi**2, pi**3, z3, z3*L2, L2, L2**2, L2**3, pi**2*L2, pi*L2**2,
         C13, C16, G, pi*G, pi*C13, pi*C16, Li3h, Ls3_60, Ls3_90, Ls3_120,
         sqrt(3)*C13, sqrt(3)*C16, log(3), pi*log(3), pi*Ls3_60]
names = ['1','pi','pi2','pi3','z3','z3log2','log2','log2^2','log2^3','pi2log2','pilog2^2',
         'Cl13','Cl16','G','piG','piCl13','piCl16','Li3h','Ls3_60','Ls3_90','Ls3_120',
         'sq3Cl13','sq3Cl16','log3','pilog3','piLs3_60']
for val, lab in [(P1,'P1'),(P1*sqrt(3),'P1*sqrt3')]:
    vv = [val]+basis
    r = pslq(vv, tol=mpf('1e-88'), maxcoeff=10**18, maxsteps=10000)
    print(lab, "rel[0]:", None if r is None else r[0])
    if r is not None and r[0]!=0:
        err = abs(sum(mpf(c)*x for c,x in zip(r, vv)))
        print("  residual:", nprint(err,4))
        out = [f"{r[0]}*T"]
        for c,n in zip(r[1:], names):
            if c: out.append(f"{c}*{n}")
        print("  "+" + ".join(out))
