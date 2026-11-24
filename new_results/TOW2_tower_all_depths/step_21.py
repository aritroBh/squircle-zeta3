import sys, os
sys.dont_write_bytecode = True
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
from towlib import *
from math import comb

def e_from_data(key, delta, s):
    i,j = key; d = i+j
    base = (2*i+1) if i==j else d
    return sum( u(i,j,base+t, delta).scale(Fr((-1)**(s-t)*comb(s,t))) for t in range(0,s+1) )

# zigzag numbers A: sec x + tan x = sum A_n x^n/n!
import sympy as sp
xx = sp.symbols('xx')
ser = sp.series(sp.sech(xx)+sp.tanh(xx), xx, 0, 16).removeO().expand()
A = {n: sp.Rational(ser.coeff(xx,n))*sp.factorial(n) for n in range(16)}
print('zigzag:', [A[n] for n in range(10)])

# check e(2s+1,s) = -i(d-s-2)*e(2s,s) for s=1..5
print('--- odd top ratio e(2s+1,s) / e(2s,s) vs -i(d-s-2)')
for s in range(1,6):
    for d in range(s+3, 12):
        key = keys_of_degree(d)[0]
        i,j = key
        base = d
        if 2*s+1 > 2*d-4 or base+s > 14: continue
        et = e_from_data(key, 2*s+1, s)
        ee = e_from_data(key, 2*s, s)
        if ee.is_zero(): continue
        r = et/ee
        conj = QSI(0,0,-(d-s-2),0)
        print(f'  s={s} d={d} key={key}: {r} vs {conj} {"OK" if r==conj else "**"}')

# alternating family F(delta,-1) = sum_s (-1)^s e(delta,s)
print('--- F(delta,-1) values')
for key in [(5,0),(6,0),(4,1)]:
    i,j=key; d=i+j
    e00 = e_from_data(key,0,0)
    print(f'key {key}:')
    for delta in range(0, 2*d-3):
        tot = ZERO
        for s in range(0, delta//2+1):
            if d+s > 14: break
            tot = tot + e_from_data(key,delta,s).scale(Fr((-1)**s))
        print(f'  d{delta}: {tot/e00 if not tot.is_zero() else ZERO}')
