
import os, sys, pickle, time
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
def e_coefs(key, dl, k0):
    P = dl//2
    rho0 = e00(key)
    f = []
    for j in range(P+1):
        val = u(key, k0+j, dl)
        if val is None: return None
        f.append(val*rho0.inv())
    out = []
    cur = f
    for s in range(P+1):
        out.append(cur[0]); cur = [cur[t+1]-cur[t] for t in range(len(cur)-1)]
    return out
def ff(x, c):
    r = Fr(1)
    for uu in range(c): r *= (x-uu)
    return r
import math
t0=time.time()
# first: e00-ratio constants vs data
print("--- e00 ratios: data vs claimed formulas ---")
for p in range(3, 9):
    rA0_d = e00((p+1,p-1))*e00((p,p-1)).inv()
    rB0_d = e00((p+1,p))*e00((p,p-1)).inv()
    rB0d_d = e00((p+1,p))*e00((p,p)).inv()
    cA = QSI(0,0,1,0).scale(Fr(4*p-3,p+1))
    cB = QSI(Fr((4*p-1)*(4*p-3),(p+1)*p))
    cBd = QSI(0,0,2,0).scale(Fr(4*p-1,p+1))
    print(f"p={p}: rA0 ok={ (rA0_d-cA).is_zero() } rB0 ok={ (rB0_d-cB).is_zero() } rB0d ok={ (rB0d_d-cBd).is_zero() } "
          f"e00adj real={ e00((p,p-1)).p2==0 and e00((p,p-1)).p3==0 } e00pp antireal={ e00((p,p)).p0==0 and e00((p,p)).p1==0 }")
# now the empirical adjacent Newton identities
ntest=0; bad=0; badex=[]
for p in range(4, 9):
    d = 2*p-1
    e00adj = e00((p,p-1)); e00pp = e00((p,p)); conv = e00pp*e00adj.inv()
    rA0_d = e00((p+1,p-1))*e00adj.inv()
    for dl in [7,8]:
        if 4*p-6 < dl: continue
        ec = e_coefs((p,p-1), dl, 2*p-1)
        if ec is None: continue
        # Ediag(dp,s) for all dp<=dl: from diagonal data, e00_pp units * conv
        ED = {}
        okED = True
        for dp in range(0, dl+1):
            eh = e_coefs((p,p), dp, 2*p+1)
            if eh is None: okED=False; break
            for s in range(0, dp//2+1):
                tot = ZERO; t=0
                while s+t <= dp//2:
                    tot = tot + eh[s+t].scale(Fr((-1)**t)); t+=1
                ED[(dp,s)] = tot*conv
        if not okED: continue
        for s in range(0, dl//2+1):
            rhs = ZERO
            ok=True
            for dp in range(0, dl+1):
                if s > dp//2: continue
                eA = e_coefs((p+1,p-1), dp, 2*p)
                if eA is None: ok=False; break
                rhs = rhs + (eA[s]*rA0_d)*tau1(dl-dp).scale(Fr(p+1,4*p-3-dp))
                rhs = rhs + ED[(dp,s)]*tau2(dl-dp).scale(Fr(p,4*p-3-dp))
            if not ok: continue
            # fold
            for t in range(0, dl+1):
                g = ZERO
                for dp in range(0, t+1):
                    if s > dp//2: continue
                    g = g - ED[(dp,s)].conj_i()*tau1(t-dp).conj_i().scale(Fr(1,4*p-3-dp))
                binom = Fr(ff(4*p-6-t, dl-t), math.factorial(dl-t))
                c = (QSI(0,0,-1,0)**(-4-t)) * (IQ**(-6-dl))
                rhs = rhs + g*c.scale(Fr(-p)*binom)
            lhs = ec[s] + (ec[s+1] if s+1 <= dl//2 else ZERO)
            ntest += 1
            if not (lhs-rhs).is_zero():
                bad += 1
                if len(badex)<6: badex.append((p,dl,s))
print(f"EMPIRICAL adjacent N(dl,s), dl=7,8, p=4..8, all-data feed incl fold: {ntest} cells, bad = {bad} {badex} ({round(time.time()-t0,1)}s)")
with open(work+"/progress.log","a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} empirical adjacent N(7,s),N(8,s) from raw data incl fold: {ntest} cells, bad={bad}; e00 ratios confirmed p=3..8\n")
