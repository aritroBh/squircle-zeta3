import json, mpmath as mp
W='<project>/research/08_theory/WT5_weight5_all_k/'
mp.mp.dps = 300
w5pp = mp.mpf(json.load(open(W+'w5pp_4800.json'))['value_full'])
E = {i: mp.mpf(json.load(open(W+f'e{i}_5800.json'))['value_full']) for i in (1,2,3)}
def binom2(x): return x*(x-1)/2
def binom3(x): return x*(x-1)*(x-2)/6
def W5poly(k):
    return w5pp + (k-5)*E[1] + binom2(k-5)*E[2] + binom3(k-5)*E[3]
# check equals previous at 5..14
from math import comb
assert all(abs(W5poly(k) - (w5pp + (k-5)*E[1] + comb(k-5,2)*E[2] + comb(k-5,3)*E[3])) < mp.mpf('1e-250') for k in range(5,15))
P4 = W5poly(4)
w5_4 = mp.mpf(json.load(open('<project>/research/08_theory/M4c_w5_closed_form/W5_5400.json'))['value_full'])
print('P(4) =', mp.nstr(P4, 52))
print('W_5(4) =', mp.nstr(w5_4, 52))
print('delta = P(4)-W_5(4) =', mp.nstr(P4-w5_4, 20))
with open(W+'progress.log','a') as f:
    f.write(f'step_M assemble: W_5(k) table k=5..14; formula vs direct quad k=6,9,14: 7.0e-210/9.0e-210/4.6e-208; k=4: delta=P(4)-W5(4)={mp.nstr(P4-w5_4,12)} != 0 (off cubic)\n')
print('logged')
