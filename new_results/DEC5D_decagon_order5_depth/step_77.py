
import pickle
import mpmath as mp
mp.mp.dps = 250
r5 = mp.sqrt(5); phi = (1+r5)/2
local250 = mp.mpf(pickle.load(open("local250.pkl","rb")))
l2 = mp.log(2); lf = mp.log(phi); z2 = mp.pi**2/6; z3 = mp.zeta(3); z5 = mp.zeta(5)
basis = {'1': mp.mpf(1), 'z5': z5, 'z2z3': z2*z3}
for a in range(6): basis[f'l2^{a}lf^{5-a}'] = l2**a*lf**(5-a)
for a in range(4): basis[f'z2 l2^{a}lf^{3-a}'] = z2*l2**a*lf**(3-a)
for a in range(3): basis[f'z3 l2^{a}lf^{2-a}'] = z3*l2**a*lf**(2-a)
basis['z2^2 l2'] = z2**2*l2; basis['z2^2 lf'] = z2**2*lf
names = list(basis)
B = [basis[n] for n in names]
vec = [local250] + B + [r5*x for x in B]
rel = mp.pslq(vec, tol=mp.mpf('1e-170'), maxcoeff=10**14, maxsteps=10000)
print("relation:", rel)
if rel:
    print("target coeff:", rel[0])
    pairs = [(names[i-1] if i>0 and i<=len(names) else '??', rel[i]) for i in range(1, len(names)+1) if rel[i]]
    qpart = [rel[len(names)+1+i] for i in range(len(names))]
    print("p-part (rational coefficients):")
    for i, n in enumerate(names):
        if rel[i+1]: print("   %-14s %s" % (n, rel[i+1]))
    print("q-part (sqrt5 coefficients):")
    for i, n in enumerate(names):
        if qpart[i]: print("   %-14s %s" % (n, qpart[i]))
