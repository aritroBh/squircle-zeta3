# step_21.py -- exact 2D tau-series through tau^12 from the derived B_2 formula; 3-way check.
import json, sympy as sp
from d16_engine import *
base='research/'
KMAX=12
B,eLam,S2 = full_S(2,KMAX)
logB2=s_log(B,KMAX)
logS2=s_add(Lambda(2,KMAX),logB2)   # log(P/(4 sqrt2))
# ---- check 1: vs H1 Route A (S_2/S_2(1) in eta=tau) ----
h1=json.load(open(base+'10_frontier/H1_octahedron_end/routeA_n2_K12.json'))['2']['eta']
zsyms={k:sp.symbols(f'z{k}') for k in range(2,15)}
def to_sympy(S,k):
    out=0
    for m,c in S.get(k,{}).items():
        t=sp.Rational(c)
        for i in m: t*=zsyms[i]
        out+=t
    return sp.expand(out)
allok=True
for k in range(13):
    diff=sp.expand(sp.sympify(h1[k]).subs({f'z{i}':zsyms[i] for i in zsyms})-to_sympy(S2,k))
    if diff!=0: allok=False; print("MISMATCH k=",k,diff)
print("[CHECK 1] S_2/S_2(1) (13 coeffs, exact, formal Z) == H1 Route A:",allok)
# ---- P/sqrt2 and log form ----
P=s_scale(S2,4)
print("\nlog(P/(4 sqrt2)) through tau^12 (single-zeta coefficients u_k(2) highlighted):")
for k in range(2,13):
    sing=[(m,c) for m,c in logS2[k].items() if len(m)==1]
    rest=[(m,c) for m,c in logS2[k].items() if len(m)>1]
    print(f" tau^{k}: single {sing} ; products {rest}")
print("\nP/sqrt(2):")
show(P,'P/sqrt2')
json.dump({'P_over_sqrt2':ser_to_json(P),'logP':ser_to_json(logS2),'B2':ser_to_json(B)}, open('tau_series_2d.json','w'), indent=1)
print("\nsaved tau_series_2d.json")
