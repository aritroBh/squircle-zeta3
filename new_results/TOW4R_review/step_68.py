
import os, sys, pickle, time
from fractions import Fraction as Fr
from math import comb
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
sys.path.insert(0, base + "/research/08_theory/TOW2_tower_all_depths")
from fastred import QSI, IQ, AQ, BQ, ZERO
H = pickle.load(open(base + "/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl",'rb'))
H.pop('_meta', None)
def q4(t): return QSI(*t)
def tau1(n):
    s = n+1
    return ((-1)**s) * ( (IQ-BQ)**(-(s+1)) + (IQ+AQ)**(-(s+1)) - IQ**(-(s+1)) )
def tau2(n):
    s = n+1
    return ((-1)**s) * ( (IQ-AQ)**(-(s+1)) + (IQ+BQ)**(-(s+1)) - IQ**(-(s+1)) )
def u(key, k, dl):
    v = H.get(key)
    return q4(v[dl][k]) if (v is not None and k in v[dl]) else None
def T1_feed(i, j, k, delta):
    d = i+j
    tot = ZERO
    for dp in range(0, delta+1):
        UA = u((i+1,j), k, dp)
        if UA is None: return None
        tot = tot + UA*tau1(delta-dp).scale(Fr(i+1, 2*d-dp-1))
        if i > j:
            UB = u((i,j+1), k, dp)
            if UB is None: return None
            tot = tot + UB*tau2(delta-dp).scale(Fr(j+1, 2*d-dp-1))
    if j == i-1:
        mm = 4*i-4-delta
        ft = ZERO
        for mu in range(mm, 4*i-3):
            g = ZERO
            for dp in range(0, 4*i-4-mu+1):
                UD = u((i,i), k, dp)
                if UD is None: return None
                g = g - UD*tau1(4*i-4-mu-dp).scale(Fr(1, 4*i-3-dp))
            ft = ft + g.conj_i() * (IQ**(-mu)) * Fr(comb(mu-2, mu-mm)) * (IQ**(mm-2))
        tot = tot + ft.scale(Fr(-i))   # -(integer i) * phi_m  (towlib: scale(-i), i = key index)
    return tot
t0=time.time()
ntest=0; bad=0; badex=[]
for d in range(12, 21):
    for p in range(d, (d-1)//2, -1):
        q = d-p
        for k in range(d+1, 53):
            for dl in [7,8]:
                if 2*d-4 < dl: continue
                got = u((p,q), k, dl)
                if got is None: continue
                pred = T1_feed(p, q, k, dl)
                if pred is None: continue
                ntest += 1
                if not (pred-got).is_zero():
                    bad += 1
                    if len(badex)<6: badex.append(((p,q),k,dl))
print(f"MY OWN T1-in-k recursion vs raw towers (FIXED fold scale -p): degrees 12..20, k=d+1..52, dl=7,8: {ntest} cells, bad = {bad} {badex} ({round(time.time()-t0,1)}s)")
with open(work+"/progress.log","a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} own T1 recursion corrected (fold scale -p integer): {ntest} cells, bad={bad}\n")
