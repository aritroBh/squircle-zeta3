
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
t0=time.time()
ntest=0; bad=0; badex=[]
for p in range(3, 14):
    e00pp = e00((p,p)); e00A = e00((p+1,p))
    conv = e00A*e00pp.inv()   # claimed rB0d = 2i(4p-1)/(p+1)
    claim = QSI(0,0,2,0).scale(Fr(4*p-1, p+1))
    assert (conv-claim).is_zero()
    for dl in [7,8]:
        if 4*p-4 < dl: continue
        eh = e_coefs((p,p), dl, 2*p+1)
        if eh is None: continue
        for s in range(0, dl//2+1):
            # claimed: ehat_pp(dl,s) = (p+1) sum_dpp tau1(dl-dpp) e_(p+1,p)(dpp,s) * rB0d /(4p-dpp-1), in e00_pp units
            tot = ZERO
            ok=True
            for dpp in range(0, dl+1):
                if s > dpp//2: continue
                eA = e_coefs((p+1,p), dpp, 2*p+1)
                if eA is None: ok=False; break
                tot = tot + (eA[s]*conv)*tau1(dl-dpp).scale(Fr(1, 4*p-dpp-1))
            if not ok: continue
            pred = tot.scale(Fr(p+1))
            ntest += 1
            if not (pred - eh[s]).is_zero():
                bad += 1
                if len(badex)<6: badex.append((p,dl,s))
print(f"EMPIRICAL diagonal ehat law, dl=7,8, p=3..13: {ntest} coefficient-cells, bad = {bad} {badex} ({round(time.time()-t0,1)}s)")
with open(work+"/progress.log","a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} empirical diagonal ehat law dl=7,8 p=3..13: {ntest} cells, bad={bad}\n")
