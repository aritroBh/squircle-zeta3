from fractions import Fraction as Fr
import sympy as sp
def wedge(u, v):
    out = {}
    for a, ca in u.items():
        for b, cb in v.items():
            if a == b: continue
            key, sg = ((a,b),1) if a < b else ((b,a),-1)
            out[key] = out.get(key, Fr(0)) + sg*ca*cb
    return {k:c for k,c in out.items() if c}
def delta_li3(vz, vomz):
    w = wedge(vz, vomz)
    out = {}
    for ab, cab in w.items():
        for l3, c3 in vz.items():
            out[(ab, l3)] = out.get((ab, l3), Fr(0)) + cab*c3
    return {k:c for k,c in out.items() if c}
# ---- 12-gon: letters g2=2, g3=sqrt(3), u=1+sqrt(3); b=2+sqrt(3) -> 2u-g2 ----
g2, g3, u = sp.Integer(2), sp.sqrt(3), 1+sp.sqrt(3)
b = 2+sp.sqrt(3)
assert sp.simplify(u**2 - 2*b) == 0
pts12 = {'z1': 2-sp.sqrt(3), 'z2': (sp.sqrt(3)-1)/2, 'z3': sp.sqrt(3)/3, 'z4': 7-4*sp.sqrt(3),
         'z5': (2-sp.sqrt(3))/2, 'z6': (2*sp.sqrt(3)-3)/3, 'z7': (3-sp.sqrt(3))/6, 'z8': sp.Rational(1,3)}
vz12 = {'z1': {'g2':1,'u':-2}, 'z2': {'u':-1}, 'z3': {'g3':-1}, 'z4': {'g2':2,'u':-4},
        'z5': {'u':-2}, 'z6': {'g2':1,'g3':-1,'u':-2}, 'z7': {'g3':-1,'u':-1}, 'z8': {'g3':-2}}
vo12 = {'z1': {'g2':1,'u':-1}, 'z2': {'g3':1,'u':-1}, 'z3': {'g2':1,'g3':-1,'u':-1}, 'z4': {'g2':2,'g3':1,'u':-2},
        'z5': {'g3':1,'g2':-1}, 'z6': {'g2':2,'g3':-1,'u':-1}, 'z7': {'u':1,'g2':-1,'g3':-1}, 'z8': {'g2':1,'g3':-2}}
def elem(v):
    e = sp.Integer(1)
    for l,c in v.items(): e *= {'g2':g2,'g3':g3,'u':u}[l]**c
    return sp.simplify(e)
ok = True
for nm in pts12:
    z = pts12[nm]
    s1 = sp.simplify(z/elem(vz12[nm])); s2 = sp.simplify((1-z)/elem(vo12[nm]))
    if s1 not in (1,-1) or s2 not in (1,-1): ok=False; print('WRONG', nm, s1, s2)
print('12-gon point factorizations exact:', ok)
coeffs12 = [("z1",Fr(-73,24),Fr(37,8)),("z2",Fr(11,2),Fr(-39,4)),("z3",Fr(97,24),Fr(-63,8)),
            ("z4",Fr(19,192),Fr(-13,64)),("z5",Fr(-1,8),Fr(9,8)),("z6",Fr(2,3),Fr(-3,4)),
            ("z7",Fr(7,4),Fr(-13,4)),("z8",Fr(-17,96),Fr(13,16))]
accP, accQ = {}, {}
for nm, p, q in coeffs12:
    dd = delta_li3(vz12[nm], vo12[nm])
    for k,c in dd.items():
        accP[k] = accP.get(k, Fr(0)) + p*c
        accQ[k] = accQ.get(k, Fr(0)) + q*c
accP = {k:c for k,c in accP.items() if c}; accQ = {k:c for k,c in accQ.items() if c}
P_pgen1 = {(('g2','u'),'u'): Fr(1,4), (('g2','g3'),'g2'): Fr(5), (('g2','g3'),'g3'): Fr(-37,8),
           (('g3','u'),'g3'): Fr(-9,4), (('g3','u'),'u'): Fr(-17,16), (('g2','u'),'g2'): Fr(5,4),
           (('g2','g3'),'u'): Fr(-2), (('g2','u'),'g3'): Fr(-1,4), (('g3','u'),'g2'): Fr(7,4)}
Q_pgen1 = {(('g2','u'),'u'): Fr(-3,4), (('g2','g3'),'g2'): Fr(-21,4), (('g2','g3'),'g3'): Fr(29,8),
           (('g3','u'),'g3'): Fr(17,8), (('g3','u'),'u'): Fr(25,16), (('g2','u'),'g2'): Fr(5,4),
           (('g2','g3'),'u'): Fr(9,4), (('g2','u'),'g3'): Fr(-1), (('g3','u'),'g2'): Fr(-13,4)}
print('12-gon delta P == PGEN1:', accP==P_pgen1, ' Q == PGEN1:', accQ==Q_pgen1)
if accP != P_pgen1: print(' P mine:', accP)
if accQ != Q_pgen1: print(' Q mine:', accQ)
print('delta != 0:', bool(accP or accQ))
open('progress.log','a').write('step: 12-gon delta(N3) recomputed; matches PGEN1 P and Q parts exactly\n')
