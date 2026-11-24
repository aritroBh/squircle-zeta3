import os, sys, pickle, json
from fractions import Fraction as Fr
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
sys.path.insert(0, workdir)
WTS = os.path.join(base, 'research/08_theory/WTS_weight_k_structure/theory')
sys.path.insert(0, WTS)
import kernel_from_tower as kft
from fastred import QSI, ZERO, pf_zero

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

# recorded kernels
wt3 = json.load(open(os.path.join(base,'research/08_theory/WT3_weight3_all_k/wt3_kernels.json')))['kernels']
def parse_kernel(rec):
    P = pf_zero()
    for pole, coeff in rec['real']:
        pr = QSI(Fr(pole[0]), Fr(pole[1]), 0, 0)
        cc = QSI(Fr(coeff[0]), Fr(coeff[1]), 0, 0)
        P[0][pr] = P[0].get(pr, ZERO) + cc
    return P

bad = []
for k in range(3, 16):
    h30 = h_orders((3,0), k)
    h21 = h_orders((2,1), k)
    K20, K11, _, _ = kft.kernel_20_21_level(h30, h21)
    rec20 = parse_kernel(wt3[str(k)]['20'])
    rec11 = parse_kernel(wt3[str(k)]['11'])
    ok20 = kft.pf_equal(K20, rec20)
    ok11 = kft.pf_equal(K11, rec11)
    if not (ok20 and ok11):
        bad.append((k, ok20, ok11))
print(f"tower->kernel vs recorded WT3 kernels k=3..15: bad = {bad}")
with open('progress.log','a') as f:
    f.write(f"step_12 tower->kernel validated vs recorded WT3 kernels k=3..15: bad={bad}\n")
