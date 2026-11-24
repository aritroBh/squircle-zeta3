
import sys, os
sys.dont_write_bytecode=True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from myengine import *
from math import comb
# Delta_pp(delta') = rho_pp(delta') - sum_s (-1)^s ehat_pp(delta',s)
def Delta(p,delta):
    out = rho(p,p,delta)
    for s in range(0,delta//2+1):
        out = out - ehat(p,delta,s).scale((-1)**s)
    return out
def Fold(p,delta):   # Fold_p(delta; Delta)
    m=4*p-4-delta; tot=ZERO
    for mu in range(m,4*p-3):
        g=ZERO
        for dp in range(0,4*p-4-mu+1):
            g = g - Delta(p,dp)*tau1(4*p-4-mu-dp).scale(Fr(1,4*p-3-dp))
        tot = tot + g.conj_i()*(IQ**(-mu))*Fr(comb(mu-2,mu-m))*(IQ**(m-2))
    return tot.scale(-p)
def feed(p,delta):
    direct=ZERO
    for dp in range(0,delta+1):
        direct = direct + Delta(p,dp)*tau2(delta-dp).scale(Fr(p,4*p-3-dp))
    return direct + Fold(p,delta)
bad=0
for p in range(2,13):
    for delta in (0,1,2):
        f = feed(p,delta)
        if not f.is_zero():
            bad+=1; print(f"E3 nonzero feed p={p} delta={delta}: {f.sympy()}")
print("E3 blindness (my encoding), p=2..12, delta=0,1,2:", "0 bad" if bad==0 else f"{bad} bad")
# also print Delta values for sanity (p=3)
for delta in (0,1,2):
    print(f"  Delta_33({delta}) =", Delta(3,delta).sympy())
