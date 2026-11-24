import sys; sys.dont_write_bytecode = True; sys.path.insert(0,'.')
from hexc2_common import cross_table, MAPS, dlog_pullback, forms_path6, FRNAM
from fractions import Fraction as Fr
tab = cross_table()
z = Fr(1,3)
e,p0,p1,pinf = tab[z][0]
print('rep (e,p0,p1,pinf):', tuple(FRNAM[q] for q in (e,p0,p1,pinf)))
e0,eh,e1,em,einf = MAPS[0]
EXP = {Fr(0):e0, Fr(1,2):eh, Fr(1):e1, Fr(-1):em, None:einf}
print('EXP:', EXP)
# correct forms: Li_1 = -int [dl(p1)-dl(pinf)]
f1 = {w:-c for w,c in dlog_pullback(EXP[p1],einf).items()}
f1 = {**f1}
from cmzv6 import add
f1 = add((f1,1), (dlog_pullback(EXP[pinf],einf), -1))
print('f1 = dl(p1)-dl(pinf):', f1)
res = forms_path6([tuple(f1.items())], EXP[p0], EXP[e])
print('int:', res, ' Li_1 = -int:', {w:-c for w,c in res.items()})
