import sys, os, pickle
sys.dont_write_bytecode=True
sys.path.insert(0, os.getcwd())
import chk_forced as C
from chk_forced import G, F, Fr, comb
H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb')); H.pop('_meta')
def QtoG(t):
    # (a,b,c,d) = (a+b r2)+i(c+d r2): only usable when b,d parts vanish -> returns (G, sqrt2 coefficient G)
    return G(t[0],t[2]), G(t[1],t[3])
def mulG(x,y): return x*y
tot=0;bad=0
for key in [(6,0),(5,1),(4,2)]:
    d=sum(key); m=key[0]
    rho=H[key][0][d]
    rho_r, rho_s = QtoG(rho)       # rho(0) = rho_r + sqrt2 * rho_s
    for dl in range(0,9):
        for k in range(d, 53):
            pred = G(0)
            for s in range(dl//2+1):
                pred = pred + F(dl,s,d,m)*comb(k-d,s)
            pr, ps = pred*rho_r, pred*rho_s       # (pred*rho): rational/sqrt2 parts
            ur, us = QtoG(H[key][dl][k])
            tot+=1
            if not (pr==ur and ps==us): bad+=1
    print(key,'checked, running bad =',bad,' cells',tot)
print('off-diagonal degree-6 towers vs own forced recursion (k=6..52, delta 0..8):',tot,'cells, bad',bad)
# negative control: corrupt tau1(8)
orig = C.tau1[8]
C.tau1[8] = G(orig.a+Fr(1,128), orig.b)
C.F.cache_clear(); C.Lam.cache_clear()
bad2=0
for (dd,mm) in [(6,6),(7,4),(9,5)]:
    if not C.closing(8,dd,mm).iszero(): bad2+=1
print('negative control (tau1(8) perturbed): closing F(8,5) fails at',bad2,'of 3 points (should be 3)')
C.tau1[8]=orig; C.F.cache_clear(); C.Lam.cache_clear()
