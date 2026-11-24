
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import math

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

# load individual reductions
reds = {}
reds.update(pickle.load(open('I_reductions_w2_7.pkl','rb')))   # w2..w7
d89 = pickle.load(open('I_reductions_w89.pkl','rb'))
reds.update(d89['reds'])                                       # w8, w9
basis8 = d89['basis8']; basis9 = d89['basis9']

# named basis element lists per weight
basis_names = {
 2: ['z2', 'log2^2'],
 3: ['z3', 'z2*log2', 'log2^3'],
 4: ['z2^2', 'z3*log2', 'z2*log2^2', 'log2^4', 'Li4'],
 5: ['z5', 'z2*z3', 'z3*log2^2', 'z2^2*log2', 'z2*log2^3', 'log2^5', 'Li5', 'Li4*log2'],
 6: ['log2^6', 'z2*log2^4', 'z2^2*log2^2', 'z2^3', 'z3*log2^3', 'z3*z2*log2', 'z3^2', 'z5*log2',
     'Li4*log2^2', 'Li4*z2', 'Li5*log2', 'Li6', 'X6'],
 7: ['log2^7', 'z2*log2^5', 'z2^2*log2^3', 'z2^3*log2', 'z3*log2^4', 'z3*z2*log2^2', 'z3*z2^2', 'z3^2*log2',
     'z5*log2^2', 'z5*z2', 'Li4*log2^3', 'Li4*z2*log2', 'Li4*z3', 'Li5*log2^2', 'Li5*z2', 'Li6*log2',
     'X6*log2', 'z7', 'Li7', 'X7a', 'X7b'],
 8: basis8,
 9: basis9,
}
# normalize names to '*' product form
def norm(nm):
    if '*' in nm: return nm
    # convert 'z2^2' to 'z2*z2', 'log2^6' to product, 'z3*log2^2' -> 'z3*log2*log2'
    parts = nm.split('*')
    out = []
    for p in parts:
        if '^' in p:
            base, e = p.split('^')
            out += [base]*int(e)
        else:
            out.append(p)
    return '*'.join(out)
basis_names = {w: [norm(x) for x in v] for w, v in basis_names.items()}

zeta_factors = {'z2','z3','z5','z7','z9'}
def is_zeta_product(nm):
    return all(f in zeta_factors for f in nm.split('*'))

e_rm = pickle.load(open('e_rm_ext.pkl','rb'))
def E(r, m): return Fr(e_rm.get(f"{r},{m}", '0'))
Dtarget = pickle.load(open('Dtilde_relations.pkl','rb'))

# map zeta-product monomials to the Dtarget basis elements (z_n single or products)
def dtarget_equiv(nm, n):
    """express zeta monomial nm in the Dtarget z-basis of weight n"""
    facs = sorted(nm.split('*'), key=lambda f: int(f[1:]))
    if n == 2: return {'z2': Fr(1)} if facs == ['z2'] else None
    if n == 3: return {'z3': Fr(1)} if facs == ['z3'] else None
    if n == 4: return {'z4': Fr(2,5)} if facs == ['z2','z2'] else None   # z4 = (2/5) z2^2
    if n == 5:
        if facs == ['z5']: return {'z5': Fr(1)}
        if facs == ['z2','z3']: return {'z2z3': Fr(1)}
    if n == 6:
        if facs == ['z2','z2','z2']: return {'z6': Fr(8,35)}
        if facs == ['z3','z3']: return {'z3^2': Fr(1)}
    if n == 7:
        if facs == ['z7']: return {'z7': Fr(1)}
        if facs == ['z2','z5']: return {'z2z5': Fr(1)}
        if facs == ['z2','z2','z3']: return {'z2^2z3': Fr(1)}
    if n == 8:
        if facs == ['z2','z2','z2','z2']: return {'z8': Fr(24,175)}
        if facs == ['z3','z5']: return {'z3z5': Fr(1)}
        if facs == ['z2','z3','z3']: return {'z2z3^2': Fr(1)}
    if n == 9:
        if facs == ['z9']: return {'z9': Fr(1)}
        if facs == ['z2','z7']: return {'z2z7': Fr(1)}
        if facs == ['z2','z2','z5']: return {'z2^2z5': Fr(1)}
        if facs == ['z2','z2','z2','z3']: return {'z2^3z3': Fr(1)}
        if facs == ['z3','z3','z3']: return {'z3^3': Fr(1)}
    return None

allok = True
for n in range(5, 10):
    bnames = basis_names[n]
    for m in range(n//2):
        acc = {}
        for j in range(1, n-2*m):
            e = E(m, n-1-j)
            if e == 0: continue
            c = 2*e/math.factorial(j)
            for cf, nm in zip(reds[(n-1-j, j)], bnames):
                if cf == 0: continue
                acc[nm] = acc.get(nm, Fr(0)) + c*cf
        # check non-zeta vanish
        nonzeta = {nm: v for nm, v in acc.items() if v and not is_zeta_product(nm)}
        # zeta part vs target
        tgt = {}
        for nm, v in acc.items():
            if v and is_zeta_product(nm):
                eq = dtarget_equiv(nm, n)
                for kk, vv in eq.items():
                    tgt[kk] = tgt.get(kk, Fr(0)) + vv*v
        tcoeffs = [Fr(c) for c in Dtarget[f"{n},{m}"]]
        keys = {'z5':'z5','z6':'z6','z7':'z7','z8':'z8','z9':'z9','z2z3':'z2z3','z2z5':'z2z5','z2z7':'z2z7',
                'z3^2':'z3^2','z3z5':'z3z5','z2z3^2':'z2z3^2','z3^3':'z3^3','z2^2z3':'z2^2z3','z2^2z5':'z2^2z5','z2^3z3':'z2^3z3'}
        # order of Dtarget basis:
        order = {5:['z5','z2z3'],6:['z6','z3^2'],7:['z7','z2z5','z2^2z3'],8:['z8','z3z5','z2z3^2'],9:['z9','z2z7','z2^2z5','z2^3z3','z3^3']}
        tgt_vec = [tgt.get(k, Fr(0)) for k in order[n]]
        match = (tgt_vec == tcoeffs)
        ok = (not nonzeta) and match
        allok &= ok
        print(f"D~_{{{n},{m}}}: non-zeta terms {len(nonzeta)}, zeta match {match} {'OK' if ok else 'PROBLEM '+str(nonzeta)+str(tgt_vec)}", flush=True)
log("assembly check (second proof route) done: %s" % ("ALL OK" if allok else "ISSUES"))
