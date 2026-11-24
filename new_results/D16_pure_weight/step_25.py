# step_24.py -- d-substitution check vs D11 a_0..a_10 (reduced); produce exact a_11, a_12; save reduced tau-series.
import json, sympy as sp
from sympy import Rational, factorial, binomial, pi, zeta
from d16_engine import *
base='research/'
KMAX=12
def from_json(J): return {int(k):{tuple(int(x) for x in m.split()):sp.Rational(cc) for m,cc in d.items()} for k,d in J.items()}
P=from_json(json.load(open('tau_series_2d.json'))['P_over_sqrt2'])
def reduce_even(expr_series,k):
    """sum monomials with even zetas reduced via sympy zeta(2j)"""
    out=0
    for m,c in expr_series.get(k,{}).items():
        t=sp.Rational(c)
        for i in m: t*=sp.zeta(i)
        out+=t
    return sp.expand(out)
# substitute tau = d/(1+d): tau^j = sum_{m>=j} (-1)^{m-j} binom(m-1,j-1) d^m
dser={}
for j in range(2,KMAX+1):
    for m in range(j,KMAX+1):
        dser[m]=dser.get(m,0)+Rational((-1)**(m-j))*binomial(m-1,j-1)*reduce_even(P,j)
dser={k:sp.expand(v) for k,v in dser.items()}
dser[0]=sp.Integer(4); dser[1]=sp.Integer(0)
txt=open(base+'08_theory/D11_rhombus_resurgence/step1_exact.txt').read()
import re
D11={}
for n in range(0,11):
    mm=re.search(rf'\na_{n} = (.*?)\n', txt, re.S)
    s=mm.group(1).strip()
    s=re.sub(r'\s+',' ',s)
    D11[n]=sp.sympify(s)
allok=True
for n in range(0,11):
    diff=sp.expand(dser.get(n,0)-D11[n])
    if diff!=0: allok=False; print("MISMATCH a_%d:"%n, diff)
print("[CHECK 3] tau-series composed with tau=d/(1+d) == D11's exact a_0..a_10:",allok)
print("\nNEW exact a_11 =", dser[11])
print("\nNEW exact a_12 =", dser[12])
print("\nnumeric a_11 =", sp.N(dser[11],30), "\nnumeric a_12 =", sp.N(dser[12],30))
# reduced tau-series (even zetas -> pi) for the report
print("\nP/sqrt(2) in tau, reduced (even zetas via pi):")
red={}
for k in range(2,13):
    red[k]=sp.expand(reduce_even(P,k))
    print(f" tau^{k}: {red[k]}")
json.dump({'a11':str(dser[11]),'a12':str(dser[12]),'tau_reduced':{str(k):str(v) for k,v in red.items()}}, open('tau_reduced_2d.json','w'), indent=1)
print("saved tau_reduced_2d.json")
