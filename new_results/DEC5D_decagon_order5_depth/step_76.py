
import pickle
import mpmath as mp
mp.mp.dps = 250
r5 = mp.sqrt(5); phi = (1+r5)/2
local250 = mp.mpf(pickle.load(open("local250.pkl","rb")))
l2 = mp.log(2); lf = mp.log(phi); z2 = mp.pi**2/6; z3 = mp.zeta(3); z5 = mp.zeta(5)
basis = {}
basis['z5'] = z5; basis['z2z3'] = z2*z3
for a in range(6):
    b = 5-a
    basis[f'l2^{a}lf^{b}'] = l2**a*lf**b
for a in range(4):
    b = 3-a
    basis[f'z2 l2^{a}lf^{b}'] = z2*l2**a*lf**b
for a in range(3):
    b = 2-a
    basis[f'z3 l2^{a}lf^{b}'] = z3*l2**a*lf**b
basis['z2^2 l2'] = z2**2*l2; basis['z2^2 lf'] = z2**2*lf
names = list(basis)
# over Q(sqrt5): represent local250/Lcoefs as p + q sqrt5 -> PSLQ on doubled system:
# find c_i in Q(sqrt5): sum c_i b_i = local. Equivalent: local = sum (p_i + q_i sqrt5) b_i.
# PSLQ on vector [local, b_1..b_n, sqrt5 b_1.. sqrt5 b_n]
B = [basis[n] for n in names]
vec = [local250] + B + [r5*x for x in B]
rel = mp.pslq(vec, tol=mp.mpf('1e-180'), maxcoeff=10**12, maxsteps=10000)
print("relation:", rel)
