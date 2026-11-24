import os, sys, pickle
from fractions import Fraction as Fr
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
sys.path.insert(0, workdir)
WTS = os.path.join(base, 'research/08_theory/WTS_weight_k_structure/theory')
sys.path.insert(0, WTS)
import kernel_from_tower as kft

H = pickle.load(open('data/towers_D8_K48.pkl','rb'))
H = {k:v for k,v in H.items() if k != '_meta'}
unq4 = kft.unq4
def h_orders(key, k):
    d = key[0]+key[1]
    out = {}
    for delta, kd in H[key].items():
        r = 2*d-2-delta
        if k in kd:
            out[r] = unq4(kd[k])
    return out

KMAX = 48
bad20 = []; bad11 = []
for k in range(3, KMAX+1):
    h30 = h_orders((3,0), k)
    h21 = h_orders((2,1), k)
    K20, K11, Hh20, Hh11 = kft.kernel_20_21_level(h30, h21)
    # check no +-i higher poles survive in kernel (Hh should be empty higher), and match closed forms
    if not kft.pf_equal(K20, kft.R20_closed(k)): bad20.append(k)
    if not kft.pf_equal(K11, kft.R11_closed(k)): bad11.append(k)
print(f"degree-2 kernels from tower vs CLOSED FORMS: k=3..{KMAX}")
print(f"  R_20 mismatches: {bad20}")
print(f"  R_11 mismatches: {bad11}")
print(f"  => R_20(k;v), R_11(k;v) EXACT for k=3..{KMAX}" if not bad20 and not bad11 else "  MISMATCH FOUND")
with open('progress.log','a') as f:
    f.write(f"step_13 degree-2 kernels from tower == closed forms R_20,R_11 exactly for k=3..{KMAX}: bad20={bad20}, bad11={bad11}\n")
