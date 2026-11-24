# B5: diagonal exceptional level k=2p and blindness E3, numerically, by MY T1 + MY candidate.
#  (i) stable extrapolation of the diagonal key to k=2p (generalized binom(-1,s)=(-1)^s) differs from the true initial rho_pp  (Delta != 0)
#  (ii) adjacent child (p,p-1) at k=2p: T1 fed with TRUE rho_pp  ==  T1 fed with STABLE extrapolation (E3 blindness), p=3..30, delta<=8
#  (iii) mutated extrapolation (shifted basis) breaks (ii)  -> test is sensitive
import sys, time
sys.dont_write_bytecode=True
from core import *
import cand
from math import comb
def ehat(p, delta, s):
    rp = rho0_tilde(p+1,p)
    eh=ZERO
    for dp in range(delta+1):
        eh = eh + TAU1[delta-dp]*cand.Fnum(dp,s,2*p+1,p+1)*fmpq(1,4*p-dp-1)
    return eh*rp*(p+1)
def stable_diag(p, delta, nshift=0):
    # sum_s ehat(delta,s) * binom(k-2p-1+nshift, s) at k=2p -> binom(-1+nshift, s) generalized
    tot=ZERO
    for s in range(delta//2+1):
        tot = tot + ehat(p,delta,s)*gbinom(-1+nshift, s)
    return tot
tot_p=0; delta_nonzero=0; e3_cells=0; e3_bad=0; mut_bad=0; mut_cells=0
for p in range(3, 31):
    d=2*p
    rho_pp = rho_tilde(p,p,8)
    ndep = min(8, 2*d-4)
    stab = [stable_diag(p,dl) for dl in range(ndep+1)]
    Dl = [rho_pp[dl]-stab[dl] for dl in range(ndep+1)]
    delta_nonzero += sum(0 if z.iszero() else 1 for z in Dl)
    tot_p+= ndep+1
    # child (p,p-1) at k = 2p : parents (p+1,p-1) [initial, degree 2p = k] and (p,p) [initial]
    rho_par = rho_tilde(p+1,p-1,8)
    def Ureal(kk,dp):
        if kk==(p+1,p-1): return rho_par[dp]
        if kk==(p,p): return rho_pp[dp]
        return ZERO
    def Ustab(kk,dp):
        if kk==(p+1,p-1): return rho_par[dp]
        if kk==(p,p): return stab[dp]
        return ZERO
    stabm=[stable_diag(p,dl,nshift=1) for dl in range(ndep+1)]
    def Umut(kk,dp):
        if kk==(p+1,p-1): return rho_par[dp]
        if kk==(p,p): return stabm[dp]
        return ZERO
    dch = 2*(2*p-1)-4
    for dl in range(min(8,dch)+1):
        a=T1_key(p,p-1,dl,Ureal); b=T1_key(p,p-1,dl,Ustab); c=T1_key(p,p-1,dl,Umut)
        e3_cells+=1
        if not (a==b): e3_bad+=1
        mut_cells+=1
        if not (a==c): mut_bad+=1
print(f'(i) diagonal (p,p), p=3..30: Delta_pp(delta)=rho_pp - stable extrapolation nonzero in {delta_nonzero}/{tot_p} cells')
print(f'(ii) E3 blindness: child (p,p-1) at k=2p, T1(true rho_pp) == T1(stable extrap): cells={e3_cells} bad={e3_bad}')
print(f'(iii) mutated extrapolation (basis shifted by one): differs from true feed in {mut_bad}/{mut_cells} cells')
