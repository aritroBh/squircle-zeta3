import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
from dod4_lib import *
d1p = pickle.load(open('dod5_task1.pkl','rb'))
TERMS = d1p['TERMS']
t0=time.time()

def Creal_acc(a, b, dps=620):
    """accurate C_real(a,b) = int_0^1 log|1-at| log|1-bt|/t dt via two methods, assert agreement."""
    mp.mp.dps = dps + 40
    pts = [mp.mpf(0), mp.mpf(1)]
    for c in (a,b):
        if c != 0:
            tc = 1/mp.mpf(c)
            if 0 < tc < 1: pts.append(tc)
    pts = sorted(pts)
    f = lambda t: mp.log(abs(1-a*t))*mp.log(abs(1-b*t))/t
    v1 = mp.quad(f, pts, method='tanh-sinh')
    v2 = mp.quad(f, pts, method='gauss-legendre', maxdegree=10)
    mp.mp.dps = dps
    d = abs(v1-v2)
    return +v1, d

worst = mp.mpf(0); res = {}
need = []
for name, kind, key, xa, ya in TERMS:
    if kind=='Lam': continue
    xf, yf = pfloat(xa), pfloat(ya)
    if not (abs(xf)<1 and abs(xf*yf)<1 and abs(yf)<=1):
        need.append((name, kind, key, xa, ya))
print(len(need), "Creal-method swaps")
for name, kind, key, xa, ya in need:
    a = mpf_of(ya); b = mpf_of(pmul(xa,ya))
    v, d = Creal_acc(a, b)
    res[(kind,key)] = (v, d)
    worst = max(worst, d)
    print("%-14s Creal=%s  |ts-gl|=%s"%(name, mp.nstr(v,15), mp.nstr(d,3)))
print("worst method-disagreement:", mp.nstr(worst,3), " t=%.0fs"%(time.time()-t0))
pickle.dump(res, open('dod5_creal_acc_tmp.pkl','wb'))
