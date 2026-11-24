from mpmath import mp, mpf, zeta, log, pi, nprint, polylog, clsin, catalan, pslq, euler
mp.dps = 160
d1 = mpf(open('d1.txt').read())
G = catalan; z3 = zeta(3); L2 = log(2); C13 = clsin(2,pi/3); C16 = clsin(2,pi/6)
Li3h = polylog(3, mpf('0.5')); Li4h = polylog(4, mpf('0.5'))
basis = [1, pi, pi**2, pi**3, zeta(2), z3, z3**2, zeta(5), L2, L2**2, L2**3,
         C13, C16, G, pi*G, pi*C13, pi*C16, Li3h, Li4h, z3*L2, pi*L2, pi**2*L2,
         pi*z3, pi**2*z3, log(3), pi*log(3), euler, pi**3*L2, zeta(3)*pi**0]
names = ['1','pi','pi2','pi3','z2','z3','z3^2','z5','log2','log2^2','log2^3','Cl13','Cl16',
         'G','piG','piCl13','piCl16','Li3h','Li4h','z3log2','pilog2','pi2log2','piz3',
         'pi2z3','log3','pilog3','gamma','pi3log2','z3b']
v = [d1] + basis
rel = pslq(v, tol=mpf('1e-145'), maxcoeff=10**22, maxsteps=10000)
print("rel[0]:", None if rel is None else rel[0])
if rel is not None and rel[0] != 0:
    err = abs(sum(mpf(c)*x for c,x in zip(rel, v)))
    print("residual:", nprint(err, 4))
    out = [f"{rel[0]}*T"]
    for c,n in zip(rel[1:], names):
        if c: out.append(f"{c}*{n}")
    print(" + ".join(out))
# also try d1/pi and d1/zeta3 quickly on smaller basis
for val,lab in [(d1/pi,'d1/pi'),(d1/z3,'d1/z3'),(d1**2,'d1^2'),(d1*pi,'d1*pi')]:
    vv = [val]+basis[:24]
    r2 = pslq(vv, tol=mpf('1e-140'), maxcoeff=10**20, maxsteps=10000)
    print(lab, "rel[0]:", None if r2 is None else r2[0], r2 if (r2 and r2[0]!=0) else '')
