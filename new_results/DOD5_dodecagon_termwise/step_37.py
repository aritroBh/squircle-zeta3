import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
from dod4_lib import *
t0=time.time()
Xs = {1: P(3,-2,3), 2: P(1,-1), 5: P(3,-1,3), 6: P(2,-1,2)}
gam={}; dl={}
IDX=[1,2,5,6]
for i in IDX:
    for j in IDX:
        if i==j: continue
        gam[(i,j)] = pdiv(Xs[j], psub(Xs[j], Xs[i]))
        dl[(i,j)]  = pmul(gam[(i,j)], psub(ONE, Xs[i]))
def C_ts(a, b, dps):
    mp.mp.dps = dps
    pts = [mp.mpf(0), mp.mpf(1)]
    for c in (a,b):
        tc = 1/mp.mpf(c)
        if 0 < tc < 1: pts.append(tc)
    pts = sorted(pts)
    return mp.quad(lambda t: mp.log(abs(1-a*t))*mp.log(abs(1-b*t))/t, pts, method='tanh-sinh')

need = {(1,2,5,'G'): (gam[(1,2)],gam[(1,5)]), (1,2,5,'D'): (dl[(1,2)],dl[(1,5)]),
        (1,2,6,'G'): (gam[(1,2)],gam[(1,6)]), (1,2,6,'D'): (dl[(1,2)],dl[(1,6)])}
out={}
for key,(pa,pb) in need.items():
    a,b = pfloat(pa), pfloat(pb)
    vs = [C_ts(a,b,dps) for dps in (620,700)]
    print(key, "C(%.7f,%.7f) = %s  |620-700|=%s  t=%.0fs"%(a,b,mp.nstr(vs[0],22),mp.nstr(abs(vs[0]-vs[1]),3),time.time()-t0))
    out[key]=vs[0]
pickle.dump(out, open('dod5_hardC_tmp.pkl','wb'))
