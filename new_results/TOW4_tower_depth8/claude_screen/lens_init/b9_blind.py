# B9: characterise the space of initial-data perturbations a_delta -> a_delta + eps*h(d,m)/Pi_{floor(delta/2)}(d) that the tower-law machinery
#     (closing + adjacent + E3, i.e. T1-evolution vs forced candidate) CANNOT detect.   Only the initial-data derivation protects against these.
import sys
sys.dont_write_bytecode=True
from core import *
import cand
eps = G(fmpq(1,1000))
def Pi(l,d):
    r=1
    for j in range(1,l+1): r*= (2*d-2*j-1)
    return r
hs = {
 '1':        lambda d,m: 1, 'd': lambda d,m: d, 'd^2': lambda d,m: d**2, 'd^3': lambda d,m: d**3, 'd^4': lambda d,m: d**4, 'd^5': lambda d,m: d**5,
 'm': lambda d,m: m, 'm^2': lambda d,m: m**2, 'd*m': lambda d,m: d*m, 'd^2*m': lambda d,m: d**2*m, 'd^3*m': lambda d,m: d**3*m,
}
def trial(delta, hname):
    h = hs[hname]; top = delta//2
    def rho_pert(p,q,N=8):
        r = list(rho_tilde(p,q,N)); d=p+q
        r[delta] = r[delta] + rho0_tilde(p,q)*eps*fmpq(h(d,p), Pi(top,d)); return r
    def base_pert(D,M):
        a=list(a_of(D,M)); a[delta]=a[delta]+eps*fmpq(h(D,M),Pi(top,D)); return a
    bad=0;n=0
    for k in range(12,22):
        cand.reset()
        T = tower_levelk(k, max(6,k-8), rho_fn=rho_pert)
        for key,vals in T.items():
            p,q=key
            if p+q==k: continue
            for dl in range(delta, 9):
                u=vals[dl]; c=cand.cand_value(p,q,dl,k,base=base_pert)
                n+=1
                if not (c==u): bad+=1
    return bad,n
for delta in (8,7):
    out=[]
    for hn in hs:
        b,n = trial(delta,hn)
        out.append(f'h={hn}: {"UNDETECTED" if b==0 else "detected ("+str(b)+"/"+str(n)+")"}')
    print(f'delta={delta}, perturbation a_delta += eps*h/Pi_{delta//2}:  ' + '; '.join(out), flush=True)
