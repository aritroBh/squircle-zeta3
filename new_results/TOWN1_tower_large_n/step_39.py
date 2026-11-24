import os, sys, pickle
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

def kernrec_to_pf(rec):
    P = pf_zero()
    for (p0,p1), c in rec['real'].items():
        P[0][QSI(p0,p1,0,0)] = QSI(*c)
    P[1] = [QSI(*c) for c in rec['ipos']]
    P[2] = [QSI(*c) for c in rec['ineg']]
    P[3] = QSI(*rec['const'])
    return P

# validate degree-3 kernels at levels k=4..18 against kern_K.pkl
bad = []
for k in range(4, 19):
    kern = pickle.load(open(os.path.join(WTS,f'kern_{k}.pkl'),'rb'))
    h40 = h_orders((4,0),k); h31 = h_orders((3,1),k); h22 = h_orders((2,2),k)
    K30, K21, _, _ = kft.kernel_deg3_level(h40, h31, h22)
    rec30 = kernrec_to_pf(kern[(3,0)])
    rec21 = kernrec_to_pf(kern[(2,1)])
    ok30 = kft.pf_equal(K30, rec30)
    ok21 = kft.pf_equal(K21, rec21)
    if not (ok30 and ok21):
        bad.append((k, ok30, ok21))
print(f"degree-3 kernels (with fold) vs recorded WTS kern_K, k=4..18: bad = {bad}")
with open('progress.log','a') as f:
    f.write(f"step_14 degree-3 kernels (fold) validated vs WTS kern_K k=4..18: bad={bad}\n")
