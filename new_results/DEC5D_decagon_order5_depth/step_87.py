
import pickle
import mpmath as mp
mp.mp.dps = 240
r5 = mp.sqrt(5); phi = (1+r5)/2
W5L = mp.mpf(pickle.load(open("w5L_250.pkl","rb"))['W5L'])
l2 = mp.log(2); lf = mp.log(phi)
z2 = mp.pi**2/6; z3 = mp.zeta(3); z5 = mp.zeta(5)
L3p = mp.polylog(3, 1/phi); L3m = mp.polylog(3, -1/phi)
L4p = mp.polylog(4, 1/phi); L4m = mp.polylog(4, -1/phi)
L5p = mp.polylog(5, 1/phi); L5m = mp.polylog(5, -1/phi)
basis = {'1': mp.mpf(1), 'z5': z5, 'z2z3': z2*z3,
         'L5p': L5p, 'L5m': L5m,
         'L4p l2': L4p*l2, 'L4p lf': L4p*lf, 'L4m l2': L4m*l2, 'L4m lf': L4m*lf,
         'L3p l2^2': L3p*l2**2, 'L3p l2lf': L3p*l2*lf, 'L3p lf^2': L3p*lf**2,
         'L3m l2^2': L3m*l2**2, 'L3m l2lf': L3m*l2*lf, 'L3m lf^2': L3m*lf**2,
         'z3 l2^2': z3*l2**2, 'z3 l2lf': z3*l2*lf, 'z3 lf^2': z3*lf**2,
         'z2^2 l2': z2**2*l2, 'z2^2 lf': z2**2*lf,
         'z2 l2^3': z2*l2**3, 'z2 l2^2lf': z2*l2**2*lf, 'z2 l2lf^2': z2*l2*lf**2, 'z2 lf^3': z2*lf**3,
         'l2^5': l2**5, 'l2^4lf': l2**4*lf, 'l2^3lf^2': l2**3*lf**2,
         'l2^2lf^3': l2**2*lf**3, 'l2lf^4': l2*lf**4, 'lf^5': lf**5}
names = list(basis)
B = [basis[n] for n in names]
vec = [W5L] + B + [r5*x for x in B]
rel = mp.pslq(vec, tol=mp.mpf('1e-160'), maxcoeff=10**12, maxsteps=20000)
print("relation found:", rel is not None)
if rel:
    print("lead:", rel[0])
    print("p-part:")
    for i, n in enumerate(names):
        if rel[i+1]: print("   %-12s %s" % (n, rel[i+1]))
    print("q-part:")
    for i, n in enumerate(names):
        j = rel[len(names)+1+i]
        if j: print("   %-12s %s" % (n, j))
    # verify
    s = sum(mp.mpf(rel[i+1])*B[i] for i in range(len(names))) + sum(mp.mpf(rel[len(names)+1+i])*r5*B[i] for i in range(len(names)))
    s = -s/mp.mpf(rel[0])
    print("residual:", mp.nstr(abs(s - W5L), 3))
