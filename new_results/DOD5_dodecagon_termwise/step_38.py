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
def C_ts(a, b, dps=660):
    mp.mp.dps = dps
    pts = [mp.mpf(0), mp.mpf(1)]
    for c in (a,b):
        tc = 1/mp.mpf(c)
        if 0 < tc < 1: pts.append(tc)
    pts = sorted(pts)
    v = mp.quad(lambda t: mp.log(abs(1-a*t))*mp.log(abs(1-b*t))/t, pts, method='tanh-sinh')
    mp.mp.dps = 620
    return +v
def SwpR_acc(a, b):
    return C_ts(a,b) + mp.log(abs(1-a))*(Li2R(b) - Li2R(b/a)) - Li3R(b)
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

vals = {}; method={}; ser_for_hard={}
for name, kind, key, xa, ya in TERMS:
    xf, yf = pfloat(xa), pfloat(ya)
    if abs(xf) < 1 and abs(xf*yf) < 1 and abs(yf) <= 1:
        vals[(kind,key)] = f_prefix_series(xa, ya); method[(kind,key)]='series'
    else:
        a = mpf_of(ya); b = mpf_of(pmul(xa,ya))
        vals[(kind,key)] = SwpR_acc(a,b); method[(kind,key)]='Creal-acc'
        if abs(xf) < 1 and abs(xf*yf) < 1:   # series converges but possibly different branch (y>1)
            ser_for_hard[(kind,key)] = f_prefix_series(xa, ya)
print("t=%.0fs, methods:"%(time.time()-t0), __import__('collections').Counter(method.values()))
# comparisons
print("\naccurate values vs stored atoms620/gdswaps620:")
worst = mp.mpf(0)
for name, kind, key, xa, ya in TERMS:
    stored = atoms['Lam'][key] if kind=='Lam' else (g['SwG'][key] if kind=='SwG' else g['SwD'][key])
    dd = abs(vals[(kind,key)] - stored)
    worst = max(worst,dd)
    if dd > mp.mpf('1e-100'):
        print("   %-14s diff=%s"%(name, mp.nstr(dd,3)))
print("worst:", mp.nstr(worst,3))
print("\nbranch check (series vs accurate real-analytic) for the y>1 convergent cases:")
for kk, sv in ser_for_hard.items():
    print("   %s: |series - realanalytic| = %s"%(kk, mp.nstr(abs(sv-vals[kk]),3)))
pickle.dump({'TERMS':TERMS,'vals':vals,'method':method,'ser_for_hard':ser_for_hard}, open('dod5_task1.pkl','wb'))
print("saved, elapsed %.0fs"%(time.time()-t0))
