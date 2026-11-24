# step_30.py -- optimal truncation point at p=1.05; Pade-Borel resummation at p=1.05 and p=1.02 vs quadrature.
import mpmath as mp, json
mp.mp.dps=130
S=[mp.mpf(x) for x in json.load(open('hi_tau_numeric_n2.json'))]
Pref={k:mp.mpf(v) for k,v in json.load(open('quad_ref.json')).items()}
from math import factorial
# Borel coefficients
Bc=[S[k]/mp.mpf(factorial(k)) for k in range(65)]
def pade_borel(tau, L, M):
    R=mp.pade(Bc[:L+M+1], L, M)
    # R is (p,q) lists; evaluate integral int_0^inf e^{-t/tau} R(t) dt
    f=lambda t: mp.exp(-t/tau)*mp.polyval(list(reversed(R[0])),t)/mp.polyval(list(reversed(R[1])),t)
    return mp.quad(f,[0,mp.inf])
for ps in ['1.02','1.05']:
    t=1-1/mp.mpf(ps); P=Pref[ps]
    print(f"--- p={ps} ---")
    # optimal truncation
    terms=[S[k]*t**k for k in range(65)]
    vals=[(abs(4*mp.sqrt(2)*sum(terms[:Nn+1])-P),Nn) for Nn in range(4,65)]
    vals.sort()
    print(f" optimal truncation: N*={vals[0][1]}, abs err={mp.nstr(vals[0][0],3)} -> {-mp.floor(mp.log10(vals[0][0]))} digits")
    for L in [20,26,30,32]:
        try:
            v=4*mp.sqrt(2)*pade_borel(t,L,L)
            print(f" Pade-Borel [{L}/{L}]: abs err={mp.nstr(abs(v-P),3)}")
        except Exception as e:
            print(f" [{L}/{L}] failed:", e)
