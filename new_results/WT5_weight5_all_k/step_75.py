import json, mpmath as mp
W='<project>/research/08_theory/WT5_weight5_all_k/'
mp.mp.dps = 300
w5pp = mp.mpf(json.load(open(W+'w5pp_4800.json'))['value_full'])
E = {i: mp.mpf(json.load(open(W+f'e{i}_5800.json'))['value_full']) for i in (1,2,3)}
from math import comb
def W5(k):
    return w5pp + (k-5)*E[1] + comb(k-5,2)*E[2] + comb(k-5,3)*E[3]
rows={}
for k in range(5,15):
    rows[k]=W5(k)
    print(k, mp.nstr(W5(k), 52))
print('--- cross-check vs direct 800-bit quadrature ---')
for k in [6,9,14]:
    q = mp.mpf(json.load(open(W+f'w5_{k}_800.json'))['value'])
    print(k, 'formula-quad =', mp.nstr(abs(W5(k)-q), 6))
print('--- k=4 off the cubic ---')
P4 = W5(4)
w5_4 = mp.mpf(json.load(open('<project>/research/08_theory/M4c_w5_closed_form/W5_5400.json'))['value_full'])
print('P(4) =', mp.nstr(P4, 52))
print('W_5(4) =', mp.nstr(w5_4, 52))
print('delta = P(4)-W_5(4) =', mp.nstr(P4-w5_4, 20))
