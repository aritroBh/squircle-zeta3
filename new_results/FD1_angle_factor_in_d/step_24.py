lib_add = r'''
# ---- squircle top table (LEAD1 result.md sec.4 closed form), beta = D/(d-1)! ----
from fractions import Fraction as Fr
from math import factorial, comb
A_Z = [1,1,1,2,5,16,61,272,1385,7936,353792,22368256,1903757312,209865342976]
def eps_m(m): return (-1)**((m+1)//2)
def delta_d(d):  # [x^d](sech+tanh) = eps_{d-1} A_d/d!
    return Fr(eps_m(d-1)*A_Z[d], factorial(d))
def squircle_D(d, i, j, pole):
    """LEAD1 closed-form top table. pole in {'-a','a','-b','0'}."""
    dd = delta_d(d)
    lastkey = (i == (d+1)//2 and j == d//2)   # the 'last key'
    diag = (i == j)
    if pole == '-a':
        return Fr((-1)**(j+1))*comb(d,j)*dd
    if pole == 'a':
        if d % 2 == 0: return Fr(0)
        v = Fr((-1)**j)*comb(d,j)*dd
        return v/2 if lastkey else v
    if pole == '-b':
        if d % 2 == 0:
            if j < d/2: return Fr((-1)**j)*comb(d,j)*dd
            return Fr(0)
        # odd d: last key also carries -b with the pole-0 value
        if lastkey: return Fr((-1)**j)*comb(d,j)*dd/2
        return Fr(0)
    if pole == '0':
        if diag or (d%2==1 and lastkey):
            return Fr((-1)**j)*comb(d,j)*dd/2
        return Fr(0)
def squircle_beta(d, i, j, pole):
    return squircle_D(d,i,j,pole)/factorial(d-1)
'''
fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
open(fd1+'/fd1_lib.py','a').write(lib_add)
import sys; sys.path.insert(0,fd1)
import importlib, fd1_lib; importlib.reload(fd1_lib)
from fd1_lib import squircle_beta, squircle_D
import pickle
bs=pickle.load(open(fd1+'/beta_squircle.pkl','rb'))
polmap={'-1':'-a','1':'a','-r':'-b','0':'0'}
bad=0; tot=0
for (d,ij,p),(v,_) in bs.items():
    tot+=1
    mine = squircle_beta(d, ij[0], ij[1], polmap[p])
    if mine != v:
        bad+=1
        if bad<=8: print("MISMATCH", d, ij, p, "table:", v, "mine:", mine)
print(f"checked {tot} squircle entries, mismatches: {bad}")
# validate d=9 against Dtop_d9.json
import json
d9=json.load(open(fd1+'/Dtop_d9.json'))
from fractions import Fraction as Fr
bad9=0
for key,poles in d9.items():
    i,j=map(int,key.split(','))
    for pole,rec in poles.items():
        ref=Fr(rec['Dtop'][0])
        mine=squircle_D(9,i,j,pole)
        if mine!=ref:
            bad9+=1; print("D9 MISMATCH", key, pole, ref, mine)
print("d=9 mismatches:", bad9)
open(fd1+'/progress.log','a').write('step_16: squircle closed-form D-table generator validated vs beta_squircle.pkl (d<=8) and Dtop_d9.json\n')
