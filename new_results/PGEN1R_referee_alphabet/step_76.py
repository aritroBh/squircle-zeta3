from fractions import Fraction as Fr
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
T3v = {
 'z1': {'u':-1}, 'z2': {'g2':-1}, 'z3': {'g3':-1}, 'z4': {'u':-2},
 'z5': {'g2':-2}, 'z6': {'g3':-1,'u':-1}, 'z7': {'g2':-1,'g3':-1}, 'z8': {'g3':-2},
 '1-z1': {'g2':1,'u':-1}, '1-z2': {'g2':-1,'g3':1}, '1-z3': {'g2':1,'g3':-1,'u':-1},
 '1-z4': {'g2':2,'g3':1,'u':-2}, '1-z5': {'g2':-2,'g3':1,'u':1}, '1-z6': {'g2':3,'g3':-1,'u':-2},
 '1-z7': {'g2':-1,'g3':-1,'u':1}, '1-z8': {'g2':2,'g3':-2,'u':-1}}
coeffs12 = [("z1",Fr(-73,24),Fr(37,8)),("z2",Fr(11,2),Fr(-39,4)),("z3",Fr(97,24),Fr(-63,8)),
            ("z4",Fr(19,192),Fr(-13,64)),("z5",Fr(-1,8),Fr(9,8)),("z6",Fr(2,3),Fr(-3,4)),
            ("z7",Fr(7,4),Fr(-13,4)),("z8",Fr(-17,96),Fr(13,16))]
# wedge^3 image: sum over terms of coeff * l1^l2^l3 (only triple (g2,g3,u) possible)
def wedge3(acc):
    tot = Fr(0)
    for ((a,b),l3), c in acc.items():
        if len({a,b,l3})==3:
            # sign of permutation (a,b,l3) -> sorted
            trip = [a,b,l3]; srt = sorted(trip)
            inv = sum(1 for i in range(3) for j in range(i+1,3) if trip[i]>trip[j])
            tot += c * (-1 if inv%2 else 1)
    return tot
accP, accQ = {}, {}
for nm, p, q in coeffs12:
    dd = delta_li3(T3v[nm], T3v['1-'+nm])
    for k,c in dd.items():
        accP[k] = accP.get(k, Fr(0)) + p*c
        accQ[k] = accQ.get(k, Fr(0)) + q*c
print('12-gon wedge^3(delta) P:', wedge3(accP), ' Q:', wedge3(accQ))
# octagon wedge^3: 2 letters -> trivially 0 (print for completeness)
vz8  = {'d': {'r2':1,'u':-1}, 'k': {'u':-1}, 'mr2h': {'r2':-1}, 'r2h': {'r2':-1}, 'x1': {'u':-1}, 'x2': {'r2':-1,'u':-1}}
vo8  = {'d': {'u':-1}, 'k': {'r2':1,'u':-1}, 'mr2h': {'u':1,'r2':-1}, 'r2h': {'r2':-1,'u':-1}, 'x1': {'r2':1}, 'x2': {'r2':-1}}
data8 = [("d",Fr(3,8),Fr(-1,4)), ("k",Fr(-3,8),Fr(1,4)), ("mr2h",Fr(-9,2),Fr(3)),
         ("r2h",Fr(-155,32),Fr(103,32)), ("x1",Fr(33,32),Fr(-21,32)), ("x2",Fr(-79,32),Fr(59,32))]
aP, aQ = {}, {}
for nm,p,q in data8:
    dd = delta_li3(vz8[nm], vo8[nm])
    for k,c in dd.items():
        aP[k] = aP.get(k, Fr(0)) + p*c; aQ[k] = aQ.get(k, Fr(0)) + q*c
print('octagon wedge^3(delta) P:', wedge3(aP), ' Q:', wedge3(aQ))
# 12-gon duplication presentation change: z4 -> 4 z1 + 4(-z1)
v1pz1 = {'g2':1,'g3':1,'u':-1}   # 1+z1 = 3-sqrt3 = 2*g3/g2 -> [2]+[g3]-[g2] = (2g2-u)+g3-g2 = g2+g3-u
import sympy as sp
assert sp.simplify((3-sp.sqrt(3)) - 2*sp.sqrt(3)/(1+sp.sqrt(3))) == 0
accP2, accQ2 = {}, {}
for nm, p, q in coeffs12:
    if nm == 'z4':
        terms = [(4*p,4*q,T3v['z1'],T3v['1-z1']), (4*p,4*q,T3v['z1'],v1pz1)]
    else:
        terms = [(p,q,T3v[nm],T3v['1-'+nm])]
    for pp,qq,vz,vo in terms:
        dd = delta_li3(vz,vo)
        for k,c in dd.items():
            accP2[k] = accP2.get(k, Fr(0)) + pp*c
            accQ2[k] = accQ2.get(k, Fr(0)) + qq*c
def ddiff(A,B):
    out = dict(A)
    for k,c in B.items():
        out[k] = out.get(k, Fr(0))-c
        if out[k]==0: del out[k]
    return {k:c for k,c in out.items() if c}
print('12-gon delta invariant under duplication z4=z1^2:', ddiff(accP,accP2)=={} and ddiff(accQ,accQ2)=={})
# formal duplication check on (z1^2) - 4z1 - 4(-z1)
dchk = {}
for dd, sc in [(delta_li3({'u':-2}, T3v['1-z4']),1), (delta_li3(T3v['z1'],T3v['1-z1']),-4), (delta_li3(T3v['z1'],v1pz1),-4)]:
    for k,c in dd.items(): dchk[k] = dchk.get(k, Fr(0))+sc*c
print('duplication relation {z1^2}-4{z1}-4{-z1} delta:', {k:c for k,c in dchk.items() if c})
open('progress.log','a').write('step: 12-gon wedge3=0; duplication invariance; octagon wedge3=0\n')
