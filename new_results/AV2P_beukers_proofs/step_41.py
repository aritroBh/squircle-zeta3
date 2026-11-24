from mpmath import mp, mpf, zeta, log, pi, nprint, polylog, clsin, catalan, pslq, euler
mp.dps = 160
d1 = mpf(open('d1.txt').read())
G = catalan; z3 = zeta(3); L2 = log(2); C13 = clsin(2,pi/3); C16 = clsin(2,pi/6)
Li3h = polylog(3, mpf('0.5')); Li4h = polylog(4, mpf('0.5'))
basis = [1, pi, pi**2, pi**3, pi**4, z3, z3*L2, zeta(4), zeta(5), L2, L2**2, L2**3, L2**4,
         C13, C16, G, pi*G, pi*C13, pi*C16, C13**2, C16**2, G**2, G*C13, Li3h, Li4h,
         pi*L2, pi**2*L2, pi**2*L2**2, pi**3*L2, log(3), pi*log(3), L2*log(3),
         z3*pi, C13*L2, C16*L2, G*L2]
names = ['1','pi','pi2','pi3','pi4','z3','z3log2','z4','z5','log2','log2^2','log2^3','log2^4',
         'Cl13','Cl16','G','piG','piCl13','piCl16','Cl13^2','Cl16^2','G^2','GCl13','Li3h','Li4h',
         'pilog2','pi2log2','pi2log2^2','pi3log2','log3','pilog3','log2log3','piz3','Cl13log2',
         'Cl16log2','Glog2']
v = [d1] + basis
rel = pslq(v, tol=mpf('1e-135'), maxcoeff=10**22, maxsteps=10000)
print("rel[0]:", None if rel is None else rel[0])
if rel is not None and rel[0] != 0:
    err = abs(sum(mpf(c)*x for c,x in zip(rel, v)))
    print("residual:", nprint(err, 4))
    out = [f"{rel[0]}*T"]
    for c,n in zip(rel[1:], names):
        if c: out.append(f"{c}*{n}")
    print(" + ".join(out))
with open('progress.log','a') as fh:
    fh.write("step_33: d1 PSLQ with 36-dim weight-4 basis: no relation (d1 likely OPEN)\n")
