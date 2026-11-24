
import sys, time, pickle, itertools
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 120
import rul1, rul1num, thx1_cert

sol9 = pickle.load(open('w9_coeffs_modp.pkl','rb'))
basis_names = pickle.load(open('I_reductions_w89.pkl','rb'))['basis9']
p0 = 4611686018427387847

def rat_recon(x, P):
    if x == 0: return Fr(0)
    r0, r1 = P, x % P
    t0, t1 = 0, 1
    while r1 * r1 * 2 > P:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        t0, t1 = t1, t0 - q * t1
    if t1 < 0: r1, t1 = -r1, -t1
    if r1 * r1 * 2 > P or t1 * t1 * 2 > P: return None
    return Fr(r1, t1)

# numeric check setup
l2 = mp.log(2)
z = {n: mp.zeta(n) for n in range(2, 11)}
Li = {n: mp.polylog(n, mp.mpf(1)/2) for n in range(4, 11)}
val = pickle.load(open('basis_vals_89.pkl','rb'))
prim_vals = dict([(nm, v) for nm, wt, v in [('log2',1,l2),('z2',2,z[2]),('z3',3,z[3]),('Li4',4,Li[4]),('z5',5,z[5]),('Li5',5,Li[5]),
    ('Li6',6,Li[6]),('X6',6,val['X6']),('z7',7,z[7]),('Li7',7,Li[7]),('X7a',7,val['X7a']),('X7b',7,val['X7b']),
    ('Li8',8,Li[8]),('X8a',8,val['X8a']),('X8b',8,val['X8b']),('X8c',8,val['X8c']),('X8d',8,val['X8d'])]])
extra_vals = {'zbar(7,1,1)': val['E9_711'], 'zbar(6,2,1)': val['E9_621'], 'zbar(5,3,1)': val['E9_531'],
              'zbar(5,1,1,1,1)': val['E9_51111'], 'zbar(4,2,1,1,1)': val['E9_42111'], 'X9fp': val['X9fp'],
              'z9': z[9], 'Li9': Li[9]}
def name_val(nm):
    if nm in extra_vals: return extra_vals[nm]
    prod = mp.mpf(1)
    for f in nm.split('*'):
        prod *= prim_vals[f]
    return prod

def E1(w, dps=110):
    return rul1num.g_word_over_0m1_num(w, dps=dps)

reds9 = {}
for (m, j), xvec in sorted(sol9.items()):
    coeffs = []
    nfail = 0
    for x in xvec:
        f = rat_recon(x, p0)
        if f is None: nfail += 1
        coeffs.append(f)
    if nfail:
        print(f"I_{{{m},{j}}}: recon fails {nfail}")
        continue
    Iv = mp.factorial(j)*mp.factorial(m)*E1((0,)*(m+1)+(-1,)*j)
    res = Iv - sum(c*name_val(nm) for c, nm in zip(coeffs, basis_names))
    reds9[(m,j)] = coeffs
    nz = sum(1 for c in coeffs if c)
    print(f"I_{{{m},{j}}}: {nz} terms, |resid| {mp.nstr(abs(res),4)}", flush=True)
pickle.dump(reds9, open('I_reductions_w9.pkl','wb'))
log("w9 reductions reconstructed + numerically checked")
