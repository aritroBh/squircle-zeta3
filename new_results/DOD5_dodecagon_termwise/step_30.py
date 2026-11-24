import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
from dod4_lib import *
t0=time.time()
mp.mp.dps = 620

IDX=[1,2,5,6]
Xs = {1: P(3,-2,3), 2: P(1,-1), 5: P(3,-1,3), 6: P(2,-1,2)}
TRIPLES = [(1,2,5),(1,2,6),(1,5,6),(2,1,5),(2,1,6),(2,5,6),(5,1,2),(5,1,6),(5,2,6),(6,1,2),(6,1,5),(6,2,5)]

gam={}; dl={}
for i in IDX:
    for j in IDX:
        if i==j: continue
        gam[(i,j)] = pdiv(Xs[j], psub(Xs[j], Xs[i]))
        dl[(i,j)]  = pmul(gam[(i,j)], psub(ONE, Xs[i]))

# ---- the 36 constants, exact arguments ----
TERMS = []   # (name, kind, key, xarg, yarg)
for i in IDX:
    for j in IDX:
        if i==j: continue
        TERMS.append(("Lam(%d,%d)"%(i,j), 'Lam', (i,j), Xs[i], pdiv(Xs[j],Xs[i])))
for (i,j,k) in TRIPLES:
    TERMS.append(("SwG(%d;%d,%d)"%(i,j,k), 'SwG', (i,j,k), pdiv(gam[(i,k)],gam[(i,j)]), gam[(i,j)]))
    TERMS.append(("SwD(%d;%d,%d)"%(i,j,k), 'SwD', (i,j,k), pdiv(dl[(i,k)],dl[(i,j)]), dl[(i,j)]))
print("total constants:", len(TERMS))

# ---- numeric evaluation ----
def f_prefix_series(x, y, dps=620):
    """f(x,y) = sum_{0<k<m} x^m y^k/(m^2 k); converges for |x|<1, |xy|<=1 (|xy|<1 unless boundary)."""
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

Li2R = lambda z: mp.re(mp.polylog(2, mp.mpc(z)))
Li3R = lambda z: mp.re(mp.polylog(3, mp.mpc(z)))
def Creal(a,b):
    return mp.quad(lambda t: mp.log(abs(1-a*t))*mp.log(abs(1-b*t))/t, [0,1])
def SwpR(a, b):
    """real-analytic f(b/a, a): C-real + ln|1-a|(Li2R(b)-Li2R(b/a)) - Li3R(b)"""
    return Creal(a,b) + mp.log(abs(1-a))*(Li2R(b) - Li2R(b/a)) - Li3R(b)
def f_value(x, y, dps=620):
    """value of the constant f(x,y): series if convergent, else real-analytic SwpR (x=b/a,y=a => a=y,b=xy)."""
    xf, yf = pfloat(x), pfloat(y)
    if abs(xf) < 1 and abs(xf*yf) < 1:
        return f_prefix_series(x, y, dps), 'series'
    mp.mp.dps = dps
    return SwpR(mpf_of(y), mpf_of(pmul(x,y))), 'Creal'

vals = {}; method = {}
worst_ser = mp.mpf(0); nser=0
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
g = pickle.load(open(D2+'gdswaps620.pkl','rb'))
for name, kind, key, xa, ya in TERMS:
    v, mth = f_value(xa, ya)
    vals[(kind,key)] = v; method[(kind,key)] = mth
print("computed all 36 values, t=%.0fs"%(time.time()-t0))
# check vs stored atoms/gdswaps
worst = mp.mpf(0)
for name, kind, key, xa, ya in TERMS:
    stored = atoms['Lam'][key] if kind=='Lam' else (g['SwG'][key] if kind=='SwG' else g['SwD'][key])
    d = abs(vals[(kind,key)] - stored)
    worst = max(worst, d)
print("all 36 fresh values vs atoms620/gdswaps620, worst:", mp.nstr(worst,3))

# series-convergence check at >=100 digits for each (count methods)
from collections import Counter
print(Counter(method.values()))
# For Creal-method swaps, additionally check the Creal formula against series in the convergent cases:
# (already series there). For the 16 Creal cases: check |x*y| etc, print table.
print("\n%-14s %-4s %-8s %s"%("name","kind","method","(x,y) approx"))
for name, kind, key, xa, ya in TERMS:
    print("%-14s %-4s %-8s x=%+.6f y=%+.6f  f=%s"%(name, kind, method[(kind,key)], pfloat(xa), pfloat(ya), mp.nstr(vals[(kind,key)],12)))
pickle.dump({'TERMS':TERMS,'vals':vals,'method':method}, open('dod5_task1.pkl','wb'))
print("saved dod5_task1.pkl, elapsed %.0fs"%(time.time()-t0))
