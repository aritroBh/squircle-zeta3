import pickle, json
from fractions import Fraction as F
from math import lcm
WORK='<project>/research/08_theory/M4f_w5_extras'
import sys
sys.dont_write_bytecode=True
sys.path.insert(0, WORK)
import nice_search as NS
sol = NS.solve(list(NS.S0), 'x')
names = sol['names']
ib2 = names.index('Li5((3-2*sqrt2))'); ibp = names.index('Li5((-1+1*sqrt2))')
ibm = names.index('Li5((1-1*sqrt2))'); im1 = names.index('Li5(-1)'); iz5 = names.index('z5')
st = pickle.loads(open(f'{WORK}/state_opt.pkl','rb').read())
bv = st['v']
v = list(bv)
for off in (0,65):
    a = v[ib2+off]
    if a: v[ibp+off] += 16*a; v[ibm+off] += 16*a; v[ib2+off] = F(0)
    a = v[im1+off]
    if a: v[iz5+off] += F(-15,16)*a; v[im1+off] = F(0)
# worst slots by numerator
items = []
slotname = lambda i: names[i] if i<65 else (names[i-65]+' [q]' if i<130 else ('Lchi' if i==130 else ('Lchi [q]' if i==131 else (f'Li5extra{(i-132)}' if i<139 else f'Li5extra{(i-139)} [q]'))))
for i,a in enumerate(v):
    if a: items.append((abs(a.numerator), a.denominator, slotname(i), a))
items.sort(reverse=True)
print('top numerator slots:')
for n,d,nm,a in items[:12]:
    print(f'  |num|={n} den={d} {nm}: {a}')
print()
# lcm contributions: factor each denominator
from sympy import factorint
D = lcm(*[a.denominator for a in v if a])
print('lcm:', D, factorint(D))
for p in [2,3,5,7]:
    best = max(( (a.denominator, slotname(i), a) for i,a in enumerate(v) if a and a.denominator % p** ( (a.denominator & -a.denominator).bit_length()-1 if p==2 else 0)==0), default=None, key=lambda t:-1)
# find slots with den divisible by 3^4 and 2^12
print('slots with den % 81 == 0 (3^4):')
for i,a in enumerate(v):
    if a and a.denominator % 81 == 0: print('   ', slotname(i), a)
print('slots with den % 4096 == 0 (2^12):')
for i,a in enumerate(v):
    if a and a.denominator % 4096 == 0: print('   ', slotname(i), a)
print('Li5/z5/Lchi coefficients:')
for i,nm in enumerate(names):
    a,b = v[i], v[65+i]
    if (a or b) and (nm.startswith('Li5') or nm=='z5'):
        print(f'  {nm:26s} {a} + ({b}) s2')
print(f'  Lchi: {v[130]} + ({v[131]}) s2')
for k in range(7):
    a,b = v[132+k], v[139+k]
    if a or b: print(f'  extra u{k}: {a} + ({b}) s2')
