
import os, sys, pickle, time, math
from fractions import Fraction as Fr
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
sys.path.insert(0, base + "/research/08_theory/TOW2_tower_all_depths")
from fastred import QSI, IQ, AQ, BQ, ZERO
H = pickle.load(open(base + "/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl",'rb'))
H.pop('_meta', None)
def q4(t): return QSI(*t)
# tau from the DL definition, exact QSI
def tau1(n):
    s = n+1
    return ((-1)**s) * ( (IQ-BQ)**(-(s+1)) + (IQ+AQ)**(-(s+1)) - IQ**(-(s+1)) )
def tau2(n):
    s = n+1
    return ((-1)**s) * ( (IQ-AQ)**(-(s+1)) + (IQ+BQ)**(-(s+1)) - IQ**(-(s+1)) )
# sanity: tau in Q(i)
for n in range(0,9):
    assert tau1(n).p1==0 and tau1(n).p3==0 and tau2(n).p1==0 and tau2(n).p3==0
def u(key, k, dl):
    v = H.get(key)
    if v is None: return None
    return q4(v[dl][k]) if k in v[dl] else None
def e00(key):
    d = sum(key); return q4(H[key][0][d])
def e_coefs(key, dl):
    """e(dl,s) = forward differences of u(dl,k)/e00 at k0=d, s=0..floor(dl/2); None if levels missing."""
    d = sum(key); k0 = d; P = dl//2
    rho0 = e00(key)
    f = []
    for j in range(P+1):
        val = u(key, k0+j, dl)
        if val is None: return None
        f.append(val*rho0.inv())
    # forward differences
    out = []
    cur = f
    for s in range(P+1):
        out.append(cur[0])
        cur = [cur[t+1]-cur[t] for t in range(len(cur)-1)]
    return out
def Lam_data(key, dl, s):
    """i(2d-1) sum_dp [tau1 eA - tau2 eB]/(2d-dp-1) from DATA coefficients."""
    p,q = key; d = p+q
    tot = ZERO
    for dp in range(0, dl+1):
        if s > dp//2: continue
        eA = e_coefs((p+1,q), dp)
        eB = e_coefs((p,q+1), dp)
        if eA is None or eB is None: return None
        t1 = tau1(dl-dp).scale(Fr(1, 2*d-dp-1))
        t2 = tau2(dl-dp).scale(Fr(1, 2*d-dp-1))
        tot = tot + eA[s]*t1 - eB[s]*t2
    return tot.scale(Fr(2*d-1))*QSI(0,0,1,0)
t0=time.time()
ntest=0; bad=0; badex=[]
for d in range(6, 21):
    for p in range(d, (d-1)//2, -1):
        q = d-p
        if p < q+2: continue   # off-diagonal only (both parents off-diagonal)
        for dl in [7,8]:
            if 2*d-4 < dl: continue
            ec = e_coefs((p,q), dl)
            if ec is None: continue
            for s in range(0, dl//2+1):
                lam = Lam_data((p,q), dl, s)
                if lam is None: continue
                lhs = ec[s] + (ec[s+1] if s+1 <= dl//2 else ZERO)
                ntest += 1
                if not (lhs-lam).is_zero():
                    bad += 1
                    if len(badex)<6: badex.append(((p,q),dl,s))
print(f"EMPIRICAL Newton N(dl,s) at dl=7,8 from raw data (off-diagonal, d=6..20): {ntest} identity-cells, bad = {bad} {badex}  ({round(time.time()-t0,1)}s)")
with open(work+"/progress.log","a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} empirical N(7,s),N(8,s) from raw TOWN1 data, d=6..20 off-diag: {ntest} cells, bad={bad}\n")
