from mpmath import mp, mpf, zeta, log, pi, nprint, polylog, clsin, catalan, pslq, quad, sqrt
mp.dps = 100
P1 = quad(lambda t: log(4/(1-t**2))**2/(3+t**2), [0,1])
G = catalan; z3 = zeta(3); L2 = log(2); C13 = clsin(2,pi/3); C16 = clsin(2,pi/6)
Li3h = polylog(3, mpf('0.5'))
Ls3_60 = quad(lambda x: log(2*mp.sin(x/2))**2, [0, pi/3])
Ls3_120 = quad(lambda x: log(2*mp.sin(x/2))**2, [0, 2*pi/3])
basis = [1, pi, pi**2, pi**3, z3, z3*L2, L2, L2**2, L2**3, pi**2*L2, pi*L2**2,
         C13, C16, G, pi*G, pi*C13, pi*C16, Li3h, Ls3_60, Ls3_120, log(3), pi*log(3),
         C13*L2, C16*L2, pi**2*C13, pi**3*L2, C13*pi**2 if False else zeta(3)*pi, sqrt(3)*C13, sqrt(3)*G]
names = ['1','pi','pi2','pi3','z3','z3log2','log2','log2^2','log2^3','pi2log2','pilog2^2',
         'Cl13','Cl16','G','piG','piCl13','piCl16','Li3h','Ls3_60','Ls3_120','log3','pilog3',
         'Cl13log2','Cl16log2','pi2Cl13','pi3log2','piz3','sq3Cl13','sq3G']
for val, lab, mm in [(P1,'P1', 10**20), (P1*sqrt(3),'P1sq3', 10**20)]:
    vv = [val]+basis
    r = pslq(vv, tol=mpf('1e-90'), maxcoeff=mm, maxsteps=10000)
    print(lab, "rel[0]:", None if r is None else r[0])
    if r is not None and r[0]!=0:
        err = abs(sum(mpf(c)*x for c,x in zip(r, vv)))
        out = [f"{r[0]}*T"]
        for c,n in zip(r[1:], names):
            if c: out.append(f"{c}*{n}")
        print("res:", nprint(err,3), " ".join(out))
with open('progress.log','a') as fh:
    fh.write("step_46: d1 = P1 - 2 zeta(3)/3 with P2=B2(1) PROVED by IBP; P1 = int log^2(4/(1-t^2))/(3+t^2) dt; PSLQ P1 open\n")
