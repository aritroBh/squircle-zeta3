from fractions import Fraction as Fr
import sympy as sp
# their letters
g2 = 1+sp.sqrt(3); g3 = sp.sqrt(3); u = 2+sp.sqrt(3)
letters12 = {'g2': g2, 'g3': g3, 'u': u}
def elem(v):
    e = sp.Integer(1)
    for l,c in v.items(): e *= letters12[l]**int(c)
    return sp.simplify(e)
T3v = {
 'z1': {'u':-1}, 'z2': {'g2':-1}, 'z3': {'g3':-1}, 'z4': {'u':-2},
 'z5': {'g2':-2}, 'z6': {'g3':-1,'u':-1}, 'z7': {'g2':-1,'g3':-1}, 'z8': {'g3':-2},
 '1-z1': {'g2':1,'u':-1}, '1-z2': {'g2':-1,'g3':1}, '1-z3': {'g2':1,'g3':-1,'u':-1},
 '1-z4': {'g2':2,'g3':1,'u':-2}, '1-z5': {'g2':-2,'g3':1,'u':1}, '1-z6': {'g2':3,'g3':-1,'u':-2},
 '1-z7': {'g2':-1,'g3':-1,'u':1}, '1-z8': {'g2':2,'g3':-2,'u':-1}}
vals = {'z1': 2-sp.sqrt(3), 'z2': (sp.sqrt(3)-1)/2, 'z3': sp.sqrt(3)/3, 'z4': 7-4*sp.sqrt(3),
        'z5': 1-sp.sqrt(3)/2, 'z6': (2*sp.sqrt(3)-3)/3, 'z7': (3-sp.sqrt(3))/6, 'z8': sp.Rational(1,3)}
ok = True
for nm, z in vals.items():
    s1 = sp.simplify(z/elem(T3v[nm])); s2 = sp.simplify((1-z)/elem(T3v['1-'+nm]))
    if s1 not in (1,-1) or s2 not in (1,-1): ok=False; print('WRONG', nm, s1, s2)
print('PGEN1 12-gon factorizations verified exact:', ok)
def wedge(uu, vv):
    out = {}
    for a, ca in uu.items():
        for bb, cb in vv.items():
            if a == bb: continue
            key, sg = ((a,bb),1) if a < bb else ((bb,a),-1)
            out[key] = out.get(key, Fr(0)) + sg*ca*cb
    return {k:c for k,c in out.items() if c}
def delta_li3(vz, vomz):
    w = wedge(vz, vomz)
    out = {}
    for ab, cab in w.items():
        for l3, c3 in vz.items():
            out[(ab, l3)] = out.get((ab, l3), Fr(0)) + cab*c3
    return {k:c for k,c in out.items() if c}
coeffs12 = [("z1",Fr(-73,24),Fr(37,8)),("z2",Fr(11,2),Fr(-39,4)),("z3",Fr(97,24),Fr(-63,8)),
            ("z4",Fr(19,192),Fr(-13,64)),("z5",Fr(-1,8),Fr(9,8)),("z6",Fr(2,3),Fr(-3,4)),
            ("z7",Fr(7,4),Fr(-13,4)),("z8",Fr(-17,96),Fr(13,16))]
accP, accQ = {}, {}
for nm, p, q in coeffs12:
    dd = delta_li3(T3v[nm], T3v['1-'+nm])
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
if accP!=P_pgen1: print('mine P:', accP)
if accQ!=Q_pgen1: print('mine Q:', accQ)
open('progress.log','a').write('step: 12-gon delta in PGEN1 letter basis; factorizations verified\n')
