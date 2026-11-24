# B7: numerical probe of the Newton/T1 identities FAR beyond the data range (exact arithmetic):
#     feed the CANDIDATE parents (forced family, my first-principles a_delta; diagonal parents: stable e-hat for k>=2p+1, TRUE initial rho_pp at k=2p)
#     into MY T1 (incl. fold) and compare with the CANDIDATE child.   Adjacent keys (p,p-1), p=3..P1; non-adjacent (d,p) p>=q+2 up to d=D1.
import sys, time
sys.dont_write_bytecode=True
from core import *
import cand
P1 = int(sys.argv[1]); D1 = int(sys.argv[2]); MODE = sys.argv[3]
_rtc={}
def _rt(p):
    if p not in _rtc: _rtc[p]=rho_tilde(p,p,8)
    return _rtc[p]
def cval(p,q,dl,k):
    if p==q and k==2*p:
        return _rt(p)[dl]
    return cand.cand_value(p,q,dl,k)
t0=time.time(); n_adj=0; bad_adj=0; n_na=0; bad_na=0; n_exc=0
for p in (range(3, P1+1) if MODE in ('adj','both') else []):
    d = 2*p-1
    for k in [d+1+j for j in range(0,10)]:      # k = 2p (exceptional diagonal level), 2p+1, ... 2p+9
        for dl in range(0, min(8, 2*d-4)+1):
            def U(kk, dp):
                return cval(kk[0], kk[1], dp, k)
            lhs = T1_key(p, p-1, dl, U)
            rhs = cand.cand_value(p, p-1, dl, k)
            n_adj+=1
            if k==2*p: n_exc+=1
            if not (lhs==rhs): bad_adj+=1
    if p%25==0: print(f'  adjacent p={p} done t={time.time()-t0:.0f}s n={n_adj} bad={bad_adj}', flush=True)
import random
random.seed(7)
for d in (range(8, D1+1) if MODE in ('na','both') else []):
    for p in range(d, d//2+1, -1):       # q = d-p, need p >= q+2  <=> p >= d/2+1
        q = d-p
        if p < q+2: continue
        if d > 60 and random.random() > 0.15: continue
        for k in [d+1+j for j in range(0,6)]:
            for dl in range(0, min(8, 2*d-4)+1):
                def U(kk, dp, _k=k):
                    return cval(kk[0], kk[1], dp, _k)
                lhs = T1_key(p, q, dl, U)
                rhs = cand.cand_value(p, q, dl, k)
                n_na+=1
                if not (lhs==rhs): bad_na+=1
    if d%20==0: print(f'  non-adjacent d={d} done t={time.time()-t0:.0f}s n={n_na} bad={bad_na}', flush=True)
print(f'B7 adjacent (p<= {P1}, d<={2*P1-1}): cells={n_adj} (of which exceptional diag level k=2p: {n_exc}) bad={bad_adj};  non-adjacent d<= {D1}: cells={n_na} bad={bad_na}; time {time.time()-t0:.0f}s')
