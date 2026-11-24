import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
from dod4_lib import *
t0=time.time()
gam12 = P(1,-1)/1  # placeholder
Xs = {1: P(3,-2,3), 2: P(1,-1), 5: P(3,-1,3), 6: P(2,-1,2)}
g12 = pdiv(Xs[2], psub(Xs[2], Xs[1]))   # gamma_{1,2}
d12 = pmul(g12, psub(ONE, Xs[1]))        # delta_{1,2}
print("gamma12 =", g12, pfloat(g12), " delta12 =", d12, pfloat(d12))

def C_ts(a, b, dps):
    mp.mp.dps = dps
    pts = [mp.mpf(0), mp.mpf(1)]
    for c in (a,b):
        tc = 1/mp.mpf(c)
        if 0 < tc < 1: pts.append(tc)
    pts = sorted(pts)
    v = mp.quad(lambda t: mp.log(abs(1-a*t))*mp.log(abs(1-b*t))/t, pts, method='tanh-sinh')
    return v

pairs = [(g12,g12),(g12,d12),(d12,d12)]
for (pa,pb) in pairs:
    a, b = pfloat(pa), pfloat(pb)
    vs = []
    for dps in (620, 680, 740):
        vs.append(C_ts(a,b,dps))
    d1 = abs(vs[0]-vs[1]); d2 = abs(vs[1]-vs[2])
    print("C(%.6f,%.6f): %s  |620-680|=%s |680-740|=%s  t=%.0fs"%(a,b,mp.nstr(vs[0],20),mp.nstr(d1,3),mp.nstr(d2,3),time.time()-t0))
