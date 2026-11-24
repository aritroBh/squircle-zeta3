import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
from dod4_lib import *
d1p = pickle.load(open('dod5_task1.pkl','rb'))
TERMS = d1p['TERMS']
t0=time.time()
mp.mp.dps = 620
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
g = pickle.load(open(D2+'gdswaps620.pkl','rb'))

Li2R = lambda z: mp.re(mp.polylog(2, mp.mpc(z)))
Li3R = lambda z: mp.re(mp.polylog(3, mp.mpc(z)))
def Creal_split(a,b):
    pts = [0,1]
    for c in (a,b):
        if c!=0 and 0 < 1/abs(c) < 1:
            # singularity where 1-c*t=0 -> t=1/c (only if in (0,1))
            tc = mp.mpf(1)/c
            if 0 < tc < 1: pts.append(tc)
    pts = sorted(set(pts))
    return mp.quad(lambda t: mp.log(abs(1-a*t))*mp.log(abs(1-b*t))/t, pts)
def SwpR2(a, b):
    return Creal_split(a,b) + mp.log(abs(1-a))*(Li2R(b) - Li2R(b/a)) - Li3R(b)
def f_prefix_series(x, y, dps=620):
    mp.mp.dps = dps+30
    x = mpf_of(x); y = mpf_of(y)
    tot = mp.mpf(0); ps = mp.mpf(0); xp = mp.mpf(1); yp = mp.mpf(1)
    for m in range(1, 2000000):
        xp *= x
        if m >= 2: tot += xp*ps/m**2
        yp *= y; ps += yp/m
        if m > 100 and abs(xp*max(1,abs(ps))) < mp.mpf(10)**(-(dps+25)): break
    mp.mp.dps = dps
    return +tot

vals = {}; method={}
for name, kind, key, xa, ya in TERMS:
    xf, yf = pfloat(xa), pfloat(ya)
    if abs(xf) < 1 and abs(xf*yf) < 1 and abs(yf) <= 1:
        vals[(kind,key)] = f_prefix_series(xa, ya); method[(kind,key)]='series'
    else:
        vals[(kind,key)] = SwpR2(mpf_of(ya), mpf_of(pmul(xa,ya))); method[(kind,key)]='Creal'
print("computed, t=%.0fs"%(time.time()-t0))
worst = mp.mpf(0); worstnm=None
for name, kind, key, xa, ya in TERMS:
    stored = atoms['Lam'][key] if kind=='Lam' else (g['SwG'][key] if kind=='SwG' else g['SwD'][key])
    dd = abs(vals[(kind,key)] - stored)
    if dd > worst: worst, worstnm = dd, name
from collections import Counter
print(Counter(method.values()))
print("all 36 vs stored, worst:", mp.nstr(worst,3), "at", worstnm)
# how many digits is that
print("digits:", -mp.log10(worst) if worst>0 else mp.inf)
pickle.dump({'TERMS':TERMS,'vals':vals,'method':method}, open('dod5_task1.pkl','wb'))
print("saved, elapsed %.0fs"%(time.time()-t0))
