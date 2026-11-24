# step_25.py -- d-substitution check vs D11 a_0..a_10 (robust parse); exact a_11,a_12; reduced tau-series.
import json, sympy as sp, re
from sympy import Rational, binomial
from d16_engine import *
base='research/'
KMAX=12
def from_json(J): return {int(k):{tuple(int(x) for x in m.split()):sp.Rational(cc) for m,cc in d.items()} for k,d in J.items()}
P=from_json(json.load(open('tau_series_2d.json'))['P_over_sqrt2'])
def red(k):
    out=0
    for m,c in P.get(k,{}).items():
        t=sp.Rational(c)
        for i in m: t*=sp.zeta(i)
        out+=t
    return sp.expand(out)
dser={}
for j in range(2,KMAX+1):
    for m in range(j,KMAX+1):
        dser[m]=dser.get(m,0)+Rational((-1)**(m-j))*binomial(m-1,j-1)*red(j)
dser={k:sp.expand(v) for k,v in dser.items()}; dser[0]=sp.Integer(4); dser[1]=sp.Integer(0)
txt=open(base+'08_theory/D11_rhombus_resurgence/step1_exact.txt').read()
D11={}
for n in range(0,11):
    i=txt.find(f'a_{n} = ')
    if i<0: print("not found",n); continue
    j=txt.find(f'a_{n+1} = ', i) if n<10 else txt.find('\n\n', i)
    if j<0: j=len(txt)
    s=' '.join(txt[i:j].split('=',1)[1].split())
    D11[n]=sp.sympify(s)
allok=True
for n in range(0,11):
    if n in D11:
        diff=sp.expand(dser.get(n,0)-D11[n])
        if diff!=0: allok=False; print("MISMATCH a_%d:"%n,diff)
print("[CHECK 3] tau-series o tau=d/(1+d) == D11 exact a_0..a_10:",allok)
print("\na_11 =", dser[11])
print("\na_12 =", dser[12])
print("\nnum a_11 =", sp.N(dser[11],25), "\nnum a_12 =", sp.N(dser[12],25))
print("\nP/sqrt2 in tau (reduced even zetas):")
redser={}
for k in range(2,13):
    redser[k]=red(k); print(f"[tau^{k}] {redser[k]}")
json.dump({'a11':str(dser[11]),'a12':str(dser[12]),'tau_reduced':{str(k):str(v) for k,v in redser.items()}}, open('tau_reduced_2d.json','w'), indent=1)
print("saved")
