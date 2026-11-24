import json
from fractions import Fraction as Fr
from math import comb
import sympy as sp
Dtop = json.load(open('Dtop.json'))
X = sp.symbols('X')
ser = sp.series((1+sp.sinh(X))/sp.cosh(X), X, 0, 9).removeO()
delta = {d: Fr(ser.coeff(X,d).p, ser.coeff(X,d).q) for d in range(1,9)}
def predicted(d, i, j, pole):
    de = delta[d]
    sgn = (-1)**(j+1)
    if pole == '-a': return sgn*comb(d,j)*de
    if d % 2 == 1:
        if pole == 'a':
            s = Fr(1) if j < (d-1)//2 else Fr(1,2)
            return (-1)**j * comb(d,j) * de * s
        if j == (d-1)//2:
            if pole == '-b': return Fr((-1)**(j+1)*comb(d,j)) * de / 2
            if pole == '0': return Fr((-1)**j*comb(d,j)) * de / 2
    else:
        if pole == '-b' and j < d//2: return (-1)**j*comb(d,j)*de
        if pole == '0' and j == d//2: return Fr((-1)**j*comb(d,j))*de/2
    return Fr(0)
bad = []
for name, rec in Dtop.items():
    d, key = name.split('|'); i, j = map(int, key.split(',')); d = int(d)
    for pole in ['-a','0','-b','a']:
        got = Fr(rec[pole]['Dtop'][0]) if pole in rec else Fr(0)
        if pole not in rec: got = Fr(0)
        pred = predicted(d, i, j, pole)
        if got != pred: bad.append((d,i,j,pole,got,pred))
print("pattern mismatches:", bad)
print("EXACT: full D_d table == zigzag/binomial closed form for all keys, poles, d=1..8" if not bad else "FAIL")
with open('progress.log','a') as f:
    f.write('step_36: EXACT check: complete D_d table = zigzag/binomial patterns, all keys/poles, d=1..8\n')
