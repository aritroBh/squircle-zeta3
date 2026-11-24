from mpmath import mp, mpf, zeta, log, pi, nprint, quad, clsin, catalan, pslq, polylog, sqrt
mp.dps = 60
J1 = quad(lambda w: log(1+w)**2/(1+w+w**2), [0,1])
J2 = quad(lambda w: log(w)*log(1+w)/(1+w+w**2), [0,1])
J3 = quad(lambda w: log(w)**2/(1+w+w**2), [0,1])
print("J1 =", nprint(J1, 30))
print("J2 =", nprint(J2, 30))
print("J3 - 8pi^3/(81 sqrt3):", nprint(J3 - 8*pi**3/(81*sqrt(3)), 4))
G = catalan; z3 = zeta(3); L2 = log(2); C13 = clsin(2,pi/3)
basis = [1, pi, pi**2, pi**3, z3, L2, L2**2, L2**3, pi**2*L2, pi*L2**2, C13, G, pi*C13, C13*L2, polylog(3,mpf('0.5')), sqrt(3)*pi**3]
names = ['1','pi','pi2','pi3','z3','log2','log2^2','log2^3','pi2log2','pilog2^2','Cl13','G','piCl13','Cl13log2','Li3h','sq3pi3']
for val, lab in [(J1,'J1'),(J2,'J2')]:
    r = pslq([val]+basis, tol=mpf('1e-52'), maxcoeff=10**14, maxsteps=10000)
    print(lab, "pslq:", r)
# check P1 assembly: P1 = (1/2)(4J1 - 4J2 + J3)
P1 = quad(lambda t: log(4/(1-t**2))**2/(3+t**2), [0,1])
print("P1 - (1/2)(4J1-4J2+J3):", nprint(P1 - (4*J1-4*J2+J3)/2, 4))
