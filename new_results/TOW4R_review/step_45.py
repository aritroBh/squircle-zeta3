
import os, sys, pickle, time, math
from fractions import Fraction as Fr
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
def e00(key): return q4(H[key][0][sum(key)])
def ff(x, c):
    r = Fr(1)
    for uu in range(c): r *= (x-uu)
    return r
t0=time.time()
ntest=0; bad=0; badex=[]
for p in range(4, 13):
    e00pp = e00((p,p))
    # Delta_pp(dp) from data for dp <= 8
    DD = {}
    ok=True
    for dp in range(0, 9):
        if 4*p-4 < dp:  # diagonal tower at depth dp needs 2(2p)-4 >= dp
            DD[dp] = None; continue
        rho_pp = u((p,p), 2*p, dp)
        eh = []
        P = dp//2
        full = rho_pp is not None
        for j in range(P+1):
            val = u((p,p), 2*p+1+j, dp)
            if val is None: full=False; break
            eh.append(val*e00pp.inv())
        if not full:
            DD[dp]=None; continue
        cur = eh; ehat=[ ]
        for s in range(P+1):
            ehat.append(cur[0]); cur=[cur[t+1]-cur[t] for t in range(len(cur)-1)]
        Dv = rho_pp*e00pp.inv()
        for s in range(P+1):
            Dv = Dv - ehat[s].scale(Fr((-1)**s))
        DD[dp] = Dv
    for dl in [7,8]:
        if 4*p-6 < dl: continue
        if any(DD[dp] is None for dp in range(0, dl+1)): continue
        feed = ZERO
        for dp in range(0, dl+1):
            feed = feed + DD[dp]*tau2(dl-dp).scale(Fr(p,4*p-3-dp))
        for t in range(0, dl+1):
            g = ZERO
            for dp in range(0, t+1):
                Dc = -DD[dp].conj_i()
                g = g - Dc*tau1(t-dp).conj_i().scale(Fr(1,4*p-3-dp))
            binom = Fr(ff(4*p-6-t, dl-t), math.factorial(dl-t))
            c = (QSI(0,0,-1,0)**(-4-t)) * (IQ**(-6-dl))
            feed = feed + g*c.scale(Fr(-p)*binom)
        ntest += 1
        if not feed.is_zero():
            bad += 1; badex.append((p,dl,feed))
print(f"EMPIRICAL E3 (diagonal blindness) dl=7,8 from raw data, p=4..12: {ntest} cells, bad = {bad} {badex[:3]} ({round(time.time()-t0,1)}s)")
with open(work+"/progress.log","a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} empirical E3 dl=7,8 from raw data p=4..12: {ntest} cells, bad={bad}\n")
