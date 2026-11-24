from mpmath import mp, mpf, zeta, log, pi, nprint, quad, polylog, clsin, catalan, pslq
mp.dps = 100
G = catalan; Cl = clsin(2, pi/3); z3 = zeta(3); L2 = log(2)
A = mpf('0.5')*quad(lambda f: polylog(2, mp.sin(f)/2)/mp.sin(f), [0, pi/2])
B = mpf('0.5')*quad(lambda f: log(mp.sin(f)/2)*log(1-mp.sin(f)/2)/mp.sin(f), [0, pi/2])
B2_2 = mpf('1.05936833579935324104250490229980343260203365771260126540342')
print("A =", nprint(A, 50))
print("B =", nprint(B, 50))
print("A+B-B2(2):", nprint(A+B-B2_2, 5))

Li3_half = polylog(3, mpf('0.5')); Li2_half = polylog(2, mpf('0.5'))
basis = [1, pi, pi**2, pi**3, z3, pi*z3, L2, pi*L2, pi**2*L2, L2**2, pi*L2**2, L2**3,
         Cl, pi*Cl, G, pi*G, Li3_half, z3*L2, Cl*L2, log(3), pi*log(3), pi**2*G,
         polylog(3, mpf('0.25')), zeta(3)*pi]
names = ['1','pi','pi2','pi3','z3','pi*z3','log2','pi*log2','pi2*log2','log2^2','pi*log2^2',
         'log2^3','Cl2(pi/3)','pi*Cl2','G','pi*G','Li3(1/2)','z3*log2','Cl2*log2','log3',
         'pi*log3','pi2*G','Li3(1/4)','z3*pi']
def hunt(target, label):
    v = [target] + basis
    rel = pslq(v, tol=mpf('1e-90'), maxcoeff=10**20, maxsteps=10000)
    print(f"--- PSLQ {label} ---")
    if rel is None: print("  no relation")
    else:
        err = abs(sum(mpf(c)*x for c,x in zip(rel, v)))
        print("  residual:", nprint(err, 4))
        out = []
        if rel[0]: out.append(f"{rel[0]}*T")
        for c,n in zip(rel[1:], names):
            if c: out.append(f"{c}*{n}")
        print("  " + " + ".join(out))
hunt(A, "A = (1/2)int Li2(sin f/2)/sin f")
hunt(B, "B = (1/2)int log log /sin f")
