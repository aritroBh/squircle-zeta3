
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import math

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

reds = {}
reds.update(pickle.load(open('I_reductions_w2_7.pkl','rb')))
d89 = pickle.load(open('I_reductions_w89.pkl','rb'))
reds.update(d89['reds'])
reds.update(pickle.load(open('I_reductions_w9.pkl','rb')))
basis8 = d89['basis8']; basis9 = d89['basis9']
basis_names = {
 5: ['z5', 'z2*z3', 'z3*log2^2', 'z2^2*log2', 'z2*log2^3', 'log2^5', 'Li5', 'Li4*log2'],
 6: ['log2^6', 'z2*log2^4', 'z2^2*log2^2', 'z2^3', 'z3*log2^3', 'z3*z2*log2', 'z3^2', 'z5*log2',
     'Li4*log2^2', 'Li4*z2', 'Li5*log2', 'Li6', 'X6'],
 7: ['log2^7', 'z2*log2^5', 'z2^2*log2^3', 'z2^3*log2', 'z3*log2^4', 'z3*z2*log2^2', 'z3*z2^2', 'z3^2*log2',
     'z5*log2^2', 'z5*z2', 'Li4*log2^3', 'Li4*z2*log2', 'Li4*z3', 'Li5*log2^2', 'Li5*z2', 'Li6*log2',
     'X6*log2', 'z7', 'Li7', 'X7a', 'X7b'],
 8: basis8,
 9: basis9,
}
def norm(nm):
    out = []
    for p in nm.split('*'):
        if '^' in p:
            b, e = p.split('^'); out += [b]*int(e)
        else: out.append(p)
    return '*'.join(out)
basis_names = {w: [norm(x) for x in v] for w, v in basis_names.items()}
zeta_factors = {'z2','z3','z5','z7','z9'}
def is_zeta_product(nm):
    return all(f in zeta_factors for f in nm.split('*'))

e_rm = pickle.load(open('e_rm_ext.pkl','rb'))
def E(r, m): return Fr(e_rm.get(f"{r},{m}", '0'))

# For each n: the I_{r,k} (r+k+1=n) reduction vectors on the named basis (columns = non-zeta parts and zeta parts).
# V_n = {c in Q^{n-1}: sum_k c_{n-1-k} red(I_{n-1-k,k}) has zero non-zeta part}
import sympy as sp
for n in range(5, 13):
    if n > 9: break
    bnames = basis_names[n]
    nonzeta_idx = [i for i, nm in enumerate(bnames) if not is_zeta_product(nm)]
    zeta_idx = [i for i, nm in enumerate(bnames) if is_zeta_product(nm)]
    # matrix: rows = I's (k=1..n-1 -> r=n-1-k), cols = non-zeta basis elements
    M = []
    for k in range(1, n):
        r = n-1-k
        cf = reds[(r, k)]
        M.append([cf[i] for i in nonzeta_idx])
    # V_n = left nullspace of M
    Ms = sp.Matrix(M)
    ns = Ms.nullspace()
    # e_m vectors
    evecs = []
    for m in range(n//2):
        ev = []
        for k in range(1, n):
            r = n-1-k
            ev.append(2*E(m, r)/math.factorial(k))
        evecs.append(sp.Matrix(ev))
    span_dim = len(ns)
    # check e_m in V_n
    ein = all((sp.Matrix(M) * ev).norm() == 0 for ev in evecs) if M else True
    # rank of e_m
    Er = sp.Matrix.hstack(*evecs).rank() if evecs else 0
    print(f"n={n}: dim V_n (zeta-landing combos) = {span_dim}, #e_m = {len(evecs)} (rank {Er}), e_m in V_n: {ein}", flush=True)
log("V_n dimensions computed from proved data")
