# step_26.py -- finalize CHECK 3 (substitute Zk->zeta(k)); inspect H1's high-order numeric data.
import json, sympy as sp
from sympy import Rational, binomial
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
base='research/'
txt=open(base+'08_theory/D11_rhombus_resurgence/step1_exact.txt').read()
D11={}
for n in range(0,11):
    i=txt.find(f'a_{n} = '); j=txt.find(f'a_{n+1} = ', i) if n<10 else txt.find('\n\n', i)
    if j<0: j=len(txt)
    s=' '.join(txt[i:j].split('=',1)[1].split())
    D11[n]=sp.sympify(s)
subs={sp.symbols(f'Z{k}'):sp.zeta(k) for k in range(2,15)}
allok=True
for n in range(0,11):
    diff=sp.expand(dser.get(n,0)-D11[n].subs(subs))
    if diff!=0: allok=False; print("MISMATCH",n,diff)
print("[CHECK 3, final] tau o (d/(1+d)) == D11 a_0..a_10 exactly:",allok)
# H1 high-order numeric file
import os
hj=json.load(open(base+'10_frontier/H1_octahedron_end/hi_coeffs_n2_K80.json'))
print("\nhi_coeffs_n2_K80.json keys:", list(hj.keys())[:5])
k0=list(hj.keys())[0]
print("sample entry:", k0, str(hj[k0])[:200])
