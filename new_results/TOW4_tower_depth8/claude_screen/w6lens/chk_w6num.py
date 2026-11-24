import sys, os, json
sys.path.insert(0, os.getcwd())
from k_load import *
import mpmath as mp
from math import comb
mp.mp.dps = 150
WT = TH + '/WT6S_weight6_symbol_all_k'
N = json.load(open(WT+'/n6_values.json'))['values']
E = json.load(open(WT+'/e5m_values.json'))['values']
print(sorted(N.keys())[:4], len(N))
POLE = {'0':(0,0),'a':(1,1),'-a':(-1,-1),'-b':(1,-1)}
r2 = mp.sqrt(2)
def W6(k):
    tot = mp.mpf(0)
    for key in KEYS5:
        for pn, pole in POLE.items():
            c = comp(k, key, pole)
            coef = mp.mpf(c[0].numerator)/c[0].denominator + r2*mp.mpf(c[1].numerator)/c[1].denominator
            tot += coef * mp.mpf(N[f'{key[0]},{key[1]}|{pn}'])
    return tot
Em = [mp.mpf(E[f'E5,{m}']) for m in range(1,6)]
worst = mp.mpf(0)
for k in range(6, 22):
    pred = sum(comb(k-6, m-1)*Em[m-1] for m in range(1,6))
    diff = abs(W6(k)-pred)
    worst = max(worst, diff)
    print(k, mp.nstr(W6(k), 25), mp.nstr(diff, 3))
print('worst |W6(k) from kernels x N  -  sum C(k-6,m-1)E5m|, k=6..21:', mp.nstr(worst,3))
# k=5 for contrast (quartic extrapolated vs actual)
def lag(ks, ys, x):
    s = mp.mpf(0)
    for i, xi in enumerate(ks):
        L = mp.mpf(1)
        for j, xj in enumerate(ks):
            if i != j: L *= mp.mpf(x - xj)/(xi - xj)
        s += ys[i]*L
    return s
ks = list(range(7, 12)); ys = [W6(k) for k in ks]
print('W6(6) actual', mp.nstr(W6(6), 30), ' quartic(7..11) at 6', mp.nstr(lag(ks, ys, 6), 30))
print('W6(5) actual', mp.nstr(W6(5), 30), ' quartic(7..11) at 5', mp.nstr(lag(ks, ys, 5), 30), ' offset', mp.nstr(W6(5)-lag(ks,ys,5), 15))
