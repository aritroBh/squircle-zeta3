import os, time, json
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
# fix delta_frac in lead2_lib.py
src = open('lead2_lib.py').read()
src = src.replace("    c = ser.coeff(x, d)\n    return Fr(int(sp.p(c)), int(sp.q(c)))",
                  "    c = sp.Rational(ser.coeff(x, d))\n    return Fr(int(c.p), int(c.q))")
open('lead2_lib.py','w').write(src)
print('patched')
import mpmath
mpmath.mp.dps = 100
from mpmath import mpf, sqrt, log, pi, quad
import sys; sys.path.insert(0,'.')
import importlib, lead2_lib
importlib.reload(lead2_lib)
from lead2_lib import load_kernels, extract_Dtop, delta_frac, predicted_D, E_coeffs, zigzag, eps
from fractions import Fraction as Fr
t0=time.time()
kern = load_kernels(18)
Dtop = extract_Dtop(kern, 9)
bad=[]
for (d,key),rec in Dtop.items():
    i,j = key
    de = delta_frac(d)
    for pole in ['-a','0','-b','a']:
        got = rec[pole]['Dtop'] if pole in rec else Fr(0)
        if got != predicted_D(d,i,j,pole,de): bad.append((d,key,pole))
poly_ok = all(pr['poly_ok'] for rec in Dtop.values() for pr in rec.values())
print("D-table d=1..9 == binomial x zigzag pattern:", not bad, " polynomial law ok:", poly_ok, f"({time.time()-t0:.0f}s)")
import pickle
pickle.dump({f"{d}|{k[0]},{k[1]}": {p: [str(pr['Dtop'].numerator), str(pr['Dtop'].denominator)] for p,pr in rec.items()} for (d,k),rec in Dtop.items()}, open('Dtop_lead2.pkl','wb'))
a = 1+sqrt(2); b = sqrt(2)-1
POLES = {'-a':-a,'0':mpf(0),'-b':-b,'a':a}
def tfun(v): return (v-b)*(v+a)/(4*v)
def ufun(v): return (a-v)*(v+b)/(4*v)
jobs = set()
for (d,key),rec in Dtop.items():
    i,j=key
    for pole,pr in rec.items():
        if pr['Dtop']!=0: jobs.add((i,j,pole))
N={}
for (i,j,pole) in sorted(jobs):
    r = POLES[pole]
    N[(i,j,pole)] = quad(lambda v: log(tfun(v))**i*log(ufun(v))**j/(v-r), [b, a])
print(f"quadrature done: {len(jobs)} integrals ({time.time()-t0:.0f}s)")
cw={}
for d in range(1,10):
    tot=mpf(0)
    for (dd,key),rec in Dtop.items():
        if dd!=d: continue
        i,j=key
        for pole,pr in rec.items():
            cf=pr['Dtop']
            if cf: tot += mpf(cf.numerator)/cf.denominator*N[(i,j,pole)]
    cw[d+1]=tot/mpmath.factorial(d-1)
E=E_coeffs(12); A=zigzag(11)
for w in range(2,11):
    m=w-2
    law=eps(m)*A[m+1]/(2*mpmath.factorial(m))*E[w]
    print(f"w={w}: |c(w)-law| = {mpmath.nstr(abs(cw[w]-law),4)}")
json.dump({str(w): mpmath.nstr(cw[w],60) for w in cw}, open('cw_dps100.json','w'))
print(f"total {time.time()-t0:.0f}s")
with open('progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: full pipeline kernels->D-table->c(w) w=2..10 dps100 vs law; saved cw_dps100.json + Dtop_lead2.pkl\n")
