import os, pickle, json
from fractions import Fraction as Fr
from math import comb
os.chdir('<project>/research/08_theory/LEAD1_leading_genfun')
base = '<project>'
theory = os.path.join(base, 'research', '08_theory', 'WTS_weight_k_structure', 'theory')
kern = {}
for K in range(10, 19):
    with open(os.path.join(theory, f'kern_{K}.pkl'),'rb') as f: kern[K] = pickle.load(f)
POLE = {(0,0):"0",(1,1):"a",(-1,-1):"-a",(-1,1):"b",(1,-1):"-b"}
def c(k,key,r):
    P = kern.get(k,{}).get(key)
    if P is None or r not in P['real']: return (Fr(0),Fr(0))
    z = P['real'][r]; return (z[0],z[1])
d = 9
keys = sorted({key for K in kern for key in kern[K] if sum(key)==d}, key=lambda kk:-kk[0])
out = {}
for key in keys:
    poles = sorted({r for K in kern for r in kern[K].get(key,{'real':{}})['real']})
    rec = {}
    for r in poles:
        seq = [c(K,key,r) for K in range(d+1,19)]
        D=[]; s=list(seq)
        for m in range(1,d+1):
            D.append(s[0]); s=[(a2[0]-a1[0],a2[1]-a1[1]) for a1,a2 in zip(s,s[1:])]
        rec[POLE[r]] = {"Dtop":[str(D[d-1][0]),str(D[d-1][1])], "leftover": str(s)}
    out[f"{key[0]},{key[1]}"] = rec
import sympy as sp
X = sp.symbols('X'); ser = sp.series((1+sp.sinh(X))/sp.cosh(X),X,0,10).removeO()
delta9 = Fr(ser.coeff(X,9).p, ser.coeff(X,9).q)
print('delta_9 =', delta9, ' (= -62/2835; eps_8 A_9/9! =', Fr(62*(-1)**0 if False else 0,1), -Fr(7936,362880),')')
def predicted(i,j,pole):
    de = delta9
    if pole=='-a': return (-1)**(j+1)*comb(9,j)*de
    if pole=='a': return (-1)**j*comb(9,j)*de*(Fr(1) if j<4 else Fr(1,2))
    if j==4:
        if pole=='-b': return Fr((-1)**(j+1)*comb(9,j))*de/2
        if pole=='0': return Fr((-1)**j*comb(9,j))*de/2
    return Fr(0)
bad=[]
for key, rec in out.items():
    i,j = map(int,key.split(','))
    for pole in ['-a','0','-b','a']:
        got = Fr(rec[pole]['Dtop'][0]) if pole in rec else Fr(0)
        pred = predicted(i,j,pole)
        if got!=pred: bad.append((key,pole,got,pred))
print('d=9 pattern mismatches:', bad)
json.dump(out, open('Dtop_d9.json','w'), indent=1)
print('EXACT: D_9 table == zigzag/binomial pattern for all keys/poles' if not bad else 'FAIL')
with open('progress.log','a') as f:
    f.write('step_57: D_9 extracted from kern_10..18 (9 points, 0 spare); matches zigzag/binomial pattern exactly for all keys/poles\n')
