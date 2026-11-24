from mpmath import mp, mpf, zeta, log, pi, nprint, polylog, clsin, catalan, pslq, quad, sqrt
mp.dps = 80
T1 = quad(lambda t: log(mp.cos(t))**2, [0, pi/6])
T2 = quad(lambda t: log(mp.sin(t))**2, [pi/6, pi/3])
T3 = quad(lambda t: log(mp.sin(t))**2, [0, pi/6])
T4 = quad(lambda t: log(mp.cos(t))*log(mp.sin(t+pi/6)), [0, pi/6])
T5 = quad(lambda t: log(mp.cos(t))*log(mp.sin(pi/6-t)), [0, pi/6])
T6 = quad(lambda t: log(mp.sin(t+pi/6))*log(mp.sin(pi/6-t)), [0, pi/6])
for n,v in [('T1',T1),('T2',T2),('T3',T3),('T4',T4),('T5',T5),('T6',T6)]:
    print(n, nprint(v, 30))
G = catalan; z3 = zeta(3); L2 = log(2); C13 = clsin(2,pi/3); C16 = clsin(2,pi/6)
Li3h = polylog(3, mpf('0.5'))
Ls3_60 = quad(lambda x: log(2*mp.sin(x/2))**2, [0, pi/3])
Ls3_30 = quad(lambda x: log(2*mp.sin(x/2))**2, [0, pi/6])
basis = [1, pi, pi**2, pi**3, z3, z3*L2, L2, L2**2, L2**3, pi**2*L2, pi*L2**2,
         C13, C16, G, pi*G, pi*C13, pi*C16, Li3h, Ls3_60, Ls3_30, log(3), pi*log(3),
         C13*L2, C16*L2, pi**2*C13]
names = ['1','pi','pi2','pi3','z3','z3log2','log2','log2^2','log2^3','pi2log2','pilog2^2',
         'Cl13','Cl16','G','piG','piCl13','piCl16','Li3h','Ls3_60','Ls3_30','log3','pilog3',
         'Cl13log2','Cl16log2','pi2Cl13']
def hunt(val, lab):
    vv = [val]+basis
    r = pslq(vv, tol=mpf('1e-72'), maxcoeff=10**16, maxsteps=10000)
    print(lab, "rel[0]:", None if r is None else r[0], end='  ')
    if r is not None and r[0]!=0:
        err = abs(sum(mpf(c)*x for c,x in zip(r, vv)))
        out = [f"{r[0]}*T"]
        for c,n in zip(r[1:], names):
            if c: out.append(f"{c}*{n}")
        print("res:", nprint(err,3), " ".join(out))
    else:
        print()
for n,v in [('T1',T1),('T2',T2),('T3',T3),('T4',T4),('T5',T5),('T6',T6)]:
    hunt(v,n)
