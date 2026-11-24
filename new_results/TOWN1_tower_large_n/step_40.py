import os, sys, pickle
from fractions import Fraction as Fr
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
sys.path.insert(0, workdir)
WTS = os.path.join(base, 'research/08_theory/WTS_weight_k_structure/theory')
sys.path.insert(0, WTS)
import kernel_from_tower as kft
from fastred import QSI, ZERO

H = pickle.load(open('data/towers_D8_K48.pkl','rb'))
H = {k:v for k,v in H.items() if k != '_meta'}
unq4 = kft.unq4
def h_orders(key, k):
    d = key[0]+key[1]
    return {2*d-2-delta: unq4(kd[k]) for delta, kd in H[key].items() if k in kd}

KMAX = 48
# collect residues per pole for K30, K21
def residues(K):
    return dict(K[0])   # {pole QSI: coeff QSI}, real poles

res30 = {}  # pole -> {k: coeff}
res21 = {}
for k in range(4, KMAX+1):
    K30, K21, Hh30, Hh21 = kft.kernel_deg3_level(h_orders((4,0),k), h_orders((3,1),k), h_orders((2,2),k))
    # confirm no +-i poles survive in kernel
    assert not [c for c in K30[1]+K30[2] if not c.is_zero()], (k,'30',K30[1],K30[2])
    assert not [c for c in K21[1]+K21[2] if not c.is_zero()], (k,'21',K21[1],K21[2])
    for pole, c in residues(K30).items():
        res30.setdefault(pole, {})[k] = c
    for pole, c in residues(K21).items():
        res21.setdefault(pole, {})[k] = c

def fdtest_deg(ys, P):
    cur = list(ys); cols=[cur]
    for j in range(1, len(ys)):
        cur=[cur[t+1]-cur[t] for t in range(len(cur)-1)]; cols.append(cur)
    le_P = all(c.is_zero() for c in cols[P+1]) if len(ys)-1>=P+1 else None
    exact = not cols[P][0].is_zero() if len(ys)-1>=P else None
    return le_P, exact

def polestr(q):
    # QSI real pole
    if q.p1==0: return f"{q.p0}"
    return f"{q.p0}{'+' if q.p1>0 else ''}{q.p1}*sqrt2"

print("K(3) kernel law test (residues polynomial in k of degree exactly 2, k=4..48):")
for name, res in [('K_30', res30), ('K_21', res21)]:
    print(f"  {name}: poles = {[polestr(p) for p in res.keys()]}")
    allok = True
    for pole, ck in sorted(res.items(), key=lambda kv: (kv[0].p0, kv[0].p1)):
        ks = sorted(ck.keys())
        ys = [ck[k] for k in ks]
        le, exact = fdtest_deg(ys, 2)
        ok = (le is True) and (exact is True)
        allok = allok and ok
        print(f"    pole {polestr(pole):>12}: deg<=2 {le}, deg==2 {exact}, pts={len(ks)}")
    print(f"  {name} ALL degree-2 polynomial: {allok}")
with open('progress.log','a') as f:
    f.write("step_15 degree-3 kernel law K(3): K_30,K_21 residues polynomial deg exactly 2 for k=4..48 (no +-i poles)\n")
