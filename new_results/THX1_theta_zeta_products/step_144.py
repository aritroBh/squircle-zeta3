
import sys, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import math
import sympy as sp

reds = {}
reds.update(pickle.load(open('I_reductions_w2_7.pkl','rb')))
d89 = pickle.load(open('I_reductions_w89.pkl','rb'))
reds.update(d89['reds'])
reds.update(pickle.load(open('I_reductions_w9.pkl','rb')))
basis_names = {
 5: ['z5', 'z2*z3', 'z3*log2^2', 'z2^2*log2', 'z2*log2^3', 'log2^5', 'Li5', 'Li4*log2'],
 7: ['log2^7', 'z2*log2^5', 'z2^2*log2^3', 'z2^3*log2', 'z3*log2^4', 'z3*z2*log2^2', 'z3*z2^2', 'z3^2*log2',
     'z5*log2^2', 'z5*z2', 'Li4*log2^3', 'Li4*z2*log2', 'Li4*z3', 'Li5*log2^2', 'Li5*z2', 'Li6*log2',
     'X6*log2', 'z7', 'Li7', 'X7a', 'X7b'],
 9: d89['basis9'],
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

for n in [5, 7, 9]:
    bnames = basis_names[n]
    nonzeta_idx = [i for i, nm in enumerate(bnames) if not is_zeta_product(nm)]
    zeta_idx = [i for i, nm in enumerate(bnames) if is_zeta_product(nm)]
    M = []
    for k in range(1, n):
        r = n-1-k
        cf = reds[(r, k)]
        M.append([cf[i] for i in nonzeta_idx])
    Ms = sp.Matrix(M)
    ns = Ms.T.nullspace()
    evecs = []
    for m in range(n//2):
        ev = []
        for k in range(1, n):
            r = n-1-k
            ev.append(2*E(m, r)/math.factorial(k))
        evecs.append(sp.Matrix(ev))
    NSe = sp.Matrix.hstack(*evecs)
    print(f"=== n={n}: V_n dim {len(ns)}")
    for v in ns:
        # is v in span(e_m)?
        aug = sp.Matrix.hstack(NSe, v)
        in_span = (aug.rank() == NSe.rank())
        # zeta part of v
        zv = {}
        for k in range(1, n):
            r = n-1-k
            cf = reds[(r, k)]
            for i in zeta_idx:
                c = v[k-1]*cf[i]
                if c:
                    zv[bnames[i]] = zv.get(bnames[i], Fr(0)) + c
        zv = {k2: x for k2, x in zv.items() if x}
        print(f"  v={list(v)} in_e_span={in_span} zeta_part={zv}")
