from mpmath import mp, mpf, zeta, log, pi, nprint, polylog, clsin, catalan, pslq, quad
mp.dps = 100
d1 = mpf(open('d1.txt').read())
# d1 = int_0^1 log^2(4/(1-t^2))/(3+t^2) dt - 2 int_0^1 log(4/(1-t^2)) t artanh t/(3+t^2) dt
P1 = quad(lambda t: log(4/(1-t**2))**2/(3+t**2), [0,1])
P2 = 2*quad(lambda t: log(4/(1-t**2))*t*mp.atanh(t)/(3+t**2), [0,1])
print("P1 =", nprint(P1, 40))
print("P2 =", nprint(P2, 40))
print("P1-P2 - d1:", nprint(P1-P2-d1, 5))
G = catalan; z3 = zeta(3); L2 = log(2); C13 = clsin(2,pi/3); C16 = clsin(2,pi/6)
Li3h = polylog(3, mpf('0.5')); Li4h = polylog(4, mpf('0.5'))
basis = [1, pi, pi**2, pi**3, pi**4, z3, z3*L2, zeta(5), L2, L2**2, L2**3, L2**4,
         C13, C16, G, pi*G, pi*C13, pi*C16, C13**2, C16**2, G**2, G*C13, Li3h, Li4h,
         pi*L2, pi**2*L2, pi**2*L2**2, pi**3*L2, log(3), pi*log(3), z3*pi, C13*L2, G*L2, pi*C13**2]
names = ['1','pi','pi2','pi3','pi4','z3','z3log2','z5','log2','log2^2','log2^3','log2^4',
         'Cl13','Cl16','G','piG','piCl13','piCl16','Cl13^2','Cl16^2','G^2','GCl13','Li3h','Li4h',
         'pilog2','pi2log2','pi2log2^2','pi3log2','log3','pilog3','piz3','Cl13log2','Glog2','piCl13^2']
for val, lab in [(P1,'P1'),(P2,'P2')]:
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
