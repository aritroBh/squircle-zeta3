import sys, os
sys.path.insert(0, os.getcwd())
from chk_feed import *
from fractions import Fraction as Fr

def lagr_vec(xs, ys, x):
    s=[Fr(0)]*4
    for i,xi in enumerate(xs):
        L=Fr(1)
        for j,xj in enumerate(xs):
            if i!=j: L*=Fr(x-xj, xi-xj)
        for t in range(4): s[t]+=ys[i][t]*L
    return tuple(s)

# discrepancy of the (3,3) tower at k=6 vs stable polynomial (fit from k=7..)
true33 = {dl: H[(3,3)][dl][6] for dl in range(9)}
stab33 = {}
for dl in range(9):
    P=dl//2; xs=list(range(7,7+P+1)); ys=[H[(3,3)][dl][k] for k in xs]
    stab33[dl]=lagr_vec(xs,ys,6)
Delta = {dl: tuple(a-b for a,b in zip(true33[dl],stab33[dl])) for dl in range(9)}

def summarize(P, label):
    simple, higher, const = pf_split(P)
    real = {p:c for p,c in realdict(simple).items() if any(x!=0 for x in c)}
    ip = [x for x in simple[1] if not x.is_zero()]
    hi = [x for x in higher[1] if not x.is_zero()]
    hi2 = [x for x in higher[2] if not x.is_zero()]
    print(f'  {label}: real-residue poles nonzero: {len(real)}; simple +-i nonzero: {len(ip)+len([x for x in simple[2] if not x.is_zero()])}; higher(+i) nonzero orders: {len(hi)}; higher(-i): {len(hi2)}; const zero: {const.is_zero()}')
    return real, ip, hi

# K32 at k=6: true tower, stable (3,3), and discrepancy-only
T_true = towers_at(6)
T_stab = towers_at(6, over33=stab33)
W_true = kernels_deg5(T_true)[(3,2)]
W_stab = kernels_deg5(T_stab)[(3,2)]
zeroT = {key:{r:QSI(0) for r in range(2,2*sum(key)-1)} for key in [(6,0),(5,1),(4,2)]}
zeroT[(3,3)] = {10-dl:QSI(*Delta[dl]) for dl in range(9)}
W_delta = kernels_deg5(zeroT)[(3,2)]
print("K32^(6): true vs stable-extrapolated (3,3) tower vs quartic")
r_true,_,_ = summarize(W_true,'true')
r_stab,_,_ = summarize(W_stab,'stable (3,3)')
print('  true real residues == stable real residues:', r_true==r_stab)
# all three keys
Wt = kernels_deg5(T_true); Ws = kernels_deg5(T_stab)
for key in KEYS5:
    a=realdict(pf_split(Wt[key])[0]); b=realdict(pf_split(Ws[key])[0])
    a={p:c for p,c in a.items() if any(c)}; b={p:c for p,c in b.items() if any(c)}
    print('  key',key,'true==stable:',a==b)
print("\nFeed of the full discrepancy Delta_33 (all depths 0..8) into W_32^(6):")
rD, ipD, hiD = summarize(W_delta,'Phi_Delta')
print('  real residues of Phi_Delta:', {(str(p[0]),str(p[1])):[str(x) for x in c] for p,c in rD.items()})
print('  higher part (+i orders) nonzero count:', len(hiD))

print("\nFeed of each single depth Delta_33(delta) into the W_32^(6) (real residues / higher part):")
for dl in range(9):
    TT = {key:{r:QSI(0) for r in range(2,2*sum(key)-1)} for key in [(6,0),(5,1),(4,2)]}
    TT[(3,3)] = {10-d2:(QSI(*Delta[d2]) if d2==dl else QSI(0)) for d2 in range(9)}
    W = kernels_deg5(TT)[(3,2)]
    simple, higher, const = pf_split(W)
    real = {p:c for p,c in realdict(simple).items() if any(x!=0 for x in c)}
    hi = [(i+1) for i,x in enumerate(higher[1]) if not x.is_zero()]
    print(f'  delta={dl}: nonzero real poles={len(real)}  higher (+i) nonzero orders={hi}  simple(+-i) nonzero={[x for x in simple[1] if not x.is_zero()]}')
