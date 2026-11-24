# B8: is the tower law (polynomial in k, degree <= floor(delta/2)) BLIND to the homogeneous perturbation a_8 -> a_8 + eps/(2d-3) ?
#     Evolve by MY T1 (with fold) from perturbed initial data (all degree-k keys incl. adjacent/diagonal) and compare with the perturbed forced candidate.
import sys
sys.dont_write_bytecode=True
from core import *
import cand
eps = G(fmpq(1,1000))
def rho_pert(p,q,N=8):
    r = list(rho_tilde(p,q,N)); d=p+q
    r[8] = r[8] + rho0_tilde(p,q)*eps/(2*d-3); return r
def base_pert(D,M):
    a=list(a_of(D,M)); a[8]=a[8]+eps/(2*D-3); return a
bad=0; n=0; badk=[]; nonadj_bad=0; adj_bad=0; diag_bad=0; n_adj=0; n_diag=0; n_off=0
for k in range(12,24):
    cand.reset()
    T = tower_levelk(k, max(6,k-8), rho_fn=rho_pert)
    for key,vals in T.items():
        p,q=key
        if p+q==k: continue
        dl=8
        u=vals[dl]
        c = cand.cand_value(p,q,dl,k,base=base_pert)
        n+=1
        isadj = (p==q+1); isdiag=(p==q)
        if isadj: n_adj+=1
        elif isdiag: n_diag+=1
        else: n_off+=1
        if not (c==u):
            bad+=1
            if isadj: adj_bad+=1
            elif isdiag: diag_bad+=1
            else: nonadj_bad+=1
            if len(badk)<6: badk.append((key,k))
print(f'delta=8 cells: {n} (off-diag non-adjacent {n_off}, adjacent {n_adj}, diagonal {n_diag}); perturbed-T1 vs perturbed-candidate mismatches: total {bad} (non-adjacent {nonadj_bad}, adjacent {adj_bad}, diagonal {diag_bad}); examples {badk}')
