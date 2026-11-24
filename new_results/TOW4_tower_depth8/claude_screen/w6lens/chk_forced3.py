import sys, os, pickle
sys.dont_write_bytecode=True
sys.path.insert(0, os.getcwd())
import chk_forced as C
from chk_forced import G, F, Fr, comb, tau1
H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb')); H.pop('_meta')
def QtoG(t): return G(t[0],t[2]), G(t[1],t[3])
p=3; d=2*p
rho_r, rho_s = QtoG(H[(3,3)][0][6])
# rho0 = rho_r + sqrt2*rho_s ; (3,3)(0) = i*sqrt2*105/16  -> rho_r = 0, rho_s = 105/16 i
rB0d = G(0,2)*(4*p-1)/(p+1)
def ehat(dl,s):
    tot=G(0)
    for dp in range(dl+1):
        tot = tot + tau1[dl-dp]*F(dp,s,2*p+1,p+1)*rB0d/(4*p-dp-1)
    return tot*(p+1)
bad=0;tot=0
for dl in range(9):
    # initial level k=6: a_delta(6,3)
    a=F(dl,0,d,p)
    ur,us=QtoG(H[(3,3)][dl][6]); tot+=1
    if not (a*rho_r==ur and a*rho_s==us): bad+=1; print('init mismatch',dl)
    for k in range(7,53):
        pred=G(0)
        for s in range(dl//2+1): pred=pred+ehat(dl,s)*comb(k-7,s)
        ur,us=QtoG(H[(3,3)][dl][k]); tot+=1
        if not (pred*rho_r==ur and pred*rho_s==us): bad+=1
print('(3,3): initial level k=6 vs a_delta(6,3) and stable levels k=7..52 vs own ehat-Newton form, delta 0..8:',tot,'cells, bad',bad)
# discrepancy Delta(dl) in rho_33(0) units, stable extrapolated to k=6: sum (-1)^s ehat
print('Delta(delta)|p=3 =',[ (lambda x: (str(x.a),str(x.b)))(F(dl,0,d,p)-sum((ehat(dl,s)*((-1)**s) for s in range(dl//2+1)),G(0))) for dl in range(9)])
