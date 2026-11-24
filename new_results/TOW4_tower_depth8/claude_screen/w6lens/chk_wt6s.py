import sys, os, json
sys.path.insert(0, os.getcwd())
from k_load import *
from math import comb
from fractions import Fraction as Fr
import sympy as sp
WT = TH + '/WT6S_weight6_symbol_all_k'
B = json.load(open(WT + '/e5m_binomial_coeffs.json'))
P = json.load(open(WT + '/deg5_residue_polynomials.json'))
print('binomial table keys:', list(B.keys()))
POLE = {'0':(0,0),'a':(1,1),'-a':(-1,-1),'-b':(1,-1)}
def diffk(key, pole, m, k0=6):
    """Delta^{m-1} of residue at k0, componentwise (4-tuple)"""
    n=m-1
    out=[Fr(0)]*4
    for j in range(n+1):
        c=(-1)**(n-j)*comb(n,j)
        v=comp(k0+j,key,pole)
        for t in range(4): out[t]+=c*v[t]
    return out
tot=0;bad=0; cells=set()
for name,ent in B.items():
    key=(int(name[0]),int(name[1])); pname=name.split('|')[1]
    pole=POLE[pname]
    assert key in KEYS5
    cells.add((key,pname))
    for m in range(1,6):
        d=diffk(key,pole,m)
        want_rat=Fr(ent['rat'][m-1]); want_s2=Fr(ent['sqrt2'][m-1])
        tot+=1
        ok = (d[0]==want_rat and d[1]==want_s2 and d[2]==0 and d[3]==0)
        if not ok:
            bad+=1; print('MISMATCH',name,m,[str(x) for x in d],'table',want_rat,want_s2)
print('(4a) difference kernels at k=6 vs e5m_binomial_coeffs.json: (key,pole,m) cells', tot, 'bad', bad, '; (key,pole) cells', len(cells))
# every (key,pole) of the kernels present in table?
miss=0
for key in KEYS5:
    for pole in POLE:
        if f'{key[0]}{key[1]}|{pole}' not in B: miss+=1; print('missing in table', key, pole)
print('table coverage missing', miss)
# the +b pole (-1,1) never present
for k in range(6,22):
    for key in KEYS5:
        assert (Fr(-1),Fr(1)) not in kern(k)[key]['real']
print('no pole at +b for k=6..21: ok')

# (4b) power-basis polynomials vs kernel for k = 6..21
kk=sp.Symbol('k')
tot=0;bad=0
for ksting,ent in P.items():
    key=(int(ksting[0]),int(ksting[1]))
    for pname,comps in ent.items():
        pole=POLE[pname]
        for k in range(6,22):
            v=comp(k,key,pole)
            rat=sp.Rational(sp.sympify(comps['rat']).subs(kk,k)) if 'rat' in comps else sp.Rational(0)
            s2=sp.Rational(sp.sympify(comps['sqrt2']).subs(kk,k)) if 'sqrt2' in comps else sp.Rational(0)
            tot+=1
            if not (Fr(int(rat.p),int(rat.q))==v[0] and Fr(int(s2.p),int(s2.q))==v[1] and v[2]==0 and v[3]==0):
                bad+=1; print('POLY MISMATCH',key,pname,k)
print('(4b) deg5_residue_polynomials.json vs kernels k=6..21: cells',tot,'bad',bad)
# (4c) polynomial degree in the json
for ksting,ent in P.items():
    for pname,comps in ent.items():
        for c,ex in comps.items():
            deg=sp.degree(sp.sympify(ex),kk)
            print(ksting,pname,c,'deg',deg, end=' | ')
    print()
