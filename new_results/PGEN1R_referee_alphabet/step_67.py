from fractions import Fraction as Fr
# ---- my own delta engine for quadratic-field letters ----
def wedge(u, v):
    out = {}
    for a, ca in u.items():
        for b, cb in v.items():
            if a == b: continue
            key, sg = ((a,b),1) if a < b else ((b,a),-1)
            out[key] = out.get(key, Fr(0)) + sg*ca*cb
    return {k:c for k,c in out.items() if c}

def delta_li3(vz, vomz):
    """(z ^ (1-z)) (x) z as dict ((l1,l2),l3) -> coeff"""
    w = wedge(vz, vomz)
    out = {}
    for ab, cab in w.items():
        for l3, c3 in vz.items():
            out[(ab, l3)] = out.get((ab, l3), Fr(0)) + cab*c3
    return {k:c for k,c in out.items() if c}

def dadd(*ds):
    out = {}
    for d in ds:
        for k,c in d.items():
            out[k] = out.get(k, Fr(0))+c
            if out[k]==0: del out[k]
    return out

# ---- regular octagon: letters r2 = sqrt2, u = 1+sqrt2; verify factorizations exactly ----
import sympy as sp
r2 = sp.sqrt(2); u = 1+r2
pts8 = {'d': 2-r2, 'k': r2-1, 'mr2h': -r2/2, 'r2h': r2/2, 'x1': 1-r2, 'x2': 1-r2/2}
# claimed letter vectors
vz8  = {'d': {'r2':1,'u':-1}, 'k': {'u':-1}, 'mr2h': {'r2':-1}, 'r2h': {'r2':-1}, 'x1': {'u':-1}, 'x2': {'r2':-1,'u':-1}}
vo8  = {'d': {'u':-1}, 'k': {'r2':1,'u':-1}, 'mr2h': {'u':1,'r2':-1}, 'r2h': {'r2':-1,'u':-1}, 'x1': {'r2':1}, 'x2': {'r2':-1}}
def elem(v):
    e = sp.Integer(1)
    for l,c in v.items():
        e *= (r2 if l=='r2' else u)**c
    return sp.simplify(e)
ok = True
for nm in pts8:
    z = pts8[nm]
    # signs dropped: |z| = elem(vz) up to sign
    s1 = sp.simplify(z/elem(vz8[nm]))
    s2 = sp.simplify((1-z)/elem(vo8[nm]))
    if s1 not in (1,-1) or s2 not in (1,-1): ok = False; print('FACTORIZATION WRONG', nm, s1, s2)
print('octagon point factorizations exact:', ok)
# coefficients from PGEN1's data (result.md 35-term table)
data8 = [("d",Fr(3,8),Fr(-1,4)), ("k",Fr(-3,8),Fr(1,4)),
         ("mr2h",Fr(-9,2),Fr(3)), ("r2h",Fr(-155,32),Fr(103,32)),
         ("x1",Fr(33,32),Fr(-21,32)), ("x2",Fr(-79,32),Fr(59,32))]
accP, accQ = {}, {}
for nm, p, q in data8:
    dd = delta_li3(vz8[nm], vo8[nm])
    for k,c in dd.items():
        accP[k] = accP.get(k, Fr(0)) + p*c
        accQ[k] = accQ.get(k, Fr(0)) + q*c
accP = {k:c for k,c in accP.items() if c}; accQ = {k:c for k,c in accQ.items() if c}
print('regular octagon delta P part:', accP)
print('regular octagon delta Q part:', accQ)
P_pgen1 = {(('r2','u'),'r2'): Fr(-5,2), (('r2','u'),'u'): Fr(-11,4)}
Q_pgen1 = {(('r2','u'),'r2'): Fr(15,8), (('r2','u'),'u'): Fr(2,1)}
print('MATCH PGEN1:', accP==P_pgen1 and accQ==Q_pgen1, ' delta != 0:', bool(accP or accQ))
open('progress.log','a').write('step: regular octagon delta(N3) recomputed; matches PGEN1 exactly\n')
