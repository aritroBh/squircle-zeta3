import os, sys, pickle, json
from fractions import Fraction as Fr
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
sys.path.insert(0, workdir)
WTS = os.path.join(base, 'research/08_theory/WTS_weight_k_structure/theory')
sys.path.insert(0, WTS)
import kernel_from_tower as kft
from fastred import QSI, ZERO, AQ, BQ

H = pickle.load(open('data/towers_D8_K48.pkl','rb'))
H = {k:v for k,v in H.items() if k != '_meta'}
unq4 = kft.unq4
def h_orders(key, k):
    d = key[0]+key[1]
    return {2*d-2-delta: unq4(kd[k]) for delta, kd in H[key].items() if k in kd}

kd3 = json.load(open(os.path.join(base,'research/08_theory/WT4_weight4_all_k/kernels_deg3.json')))
polesym = {'0': ZERO, 'a': AQ, '-a': -AQ, 'b': BQ, '-b': -BQ}
def rec_to_pf(rec):
    P = { }
    for pole, coeff in rec:
        P[polesym[pole]] = QSI(Fr(coeff[0]), Fr(coeff[1]), 0, 0)
    return P

bad=[]
for k in range(4, 13):
    K30, K21, _, _ = kft.kernel_deg3_level(h_orders((4,0),k), h_orders((3,1),k), h_orders((2,2),k))
    r30 = rec_to_pf(kd3[str(k)]['3,0'])
    r21 = rec_to_pf(kd3[str(k)]['2,1'])
    ok30 = all(K30[0].get(p,ZERO)==c for p,c in r30.items()) and all(c.is_zero() for p,c in K30[0].items() if p not in r30)
    ok21 = all(K21[0].get(p,ZERO)==c for p,c in r21.items()) and all(c.is_zero() for p,c in K21[0].items() if p not in r21)
    if not (ok30 and ok21): bad.append((k,ok30,ok21))
print("degree-3 kernels vs WT4 kernels_deg3.json k=4..12: bad =", bad)
with open('progress.log','a') as f:
    f.write(f"step_16 degree-3 kernels match WT4 kernels_deg3.json k=4..12: bad={bad}\n")

# Extract explicit degree-2 polynomial coefficients for each residue (rational & sqrt2 parts separately)
def fit_poly2(ks, ys):
    # returns (c0,c1,c2) rational with y = c0 + c1*k + c2*k^2, using first 3 pts, exact
    from math import comb
    cur=list(ys); cols=[cur]
    for j in range(1,3):
        cur=[cur[t+1]-cur[t] for t in range(len(cur)-1)]; cols.append(cur)
    k0=ks[0]
    # y(k) = cols[0][0] + cols[1][0]*(k-k0) + cols[2][0]/2*(k-k0)(k-k0-1)
    d0,d1,d2 = cols[0][0], cols[1][0], cols[2][0]/2
    # expand in k
    # (k-k0) = k-k0 ; (k-k0)(k-k0-1) = k^2 -(2k0+1)k + k0(k0+1)
    c2 = d2
    c1 = d1 - d2*Fr(2*k0+1)
    c0 = d0 - d1*k0 + d2*Fr(k0*(k0+1))
    return (c0,c1,c2)

def res_polys(keyname):
    # collect per pole the (rational, sqrt2) coefficient sequences k=4..48
    seqs = {}
    for k in range(4,49):
        K30, K21, _, _ = kft.kernel_deg3_level(h_orders((4,0),k), h_orders((3,1),k), h_orders((2,2),k))
        K = K30 if keyname=='30' else K21
        for p,c in K[0].items():
            seqs.setdefault(p,{})[k]=c
    out={}
    for p,ck in seqs.items():
        ks=sorted(ck)
        ys=[ck[k] for k in ks]
        rat=[c.scale(1) for c in ys]  # full QSI but real: p2=p3=0
        # fit on the QSI (real) values: use p0,p1 separately
        p0=[c.p0 for c in ys]; p1=[c.p1 for c in ys]
        out[p]=(fit_poly2(ks,p0), fit_poly2(ks,p1))
    return out

for nm in ['30','21']:
    print(f"K_{nm} residue polynomials (y = c0 + c1 k + c2 k^2), exact for k=4..48:")
    for p,(q0,q1) in res_polys(nm).items():
        ps = {'0':'0','a':'a','-a':'-a','b':'b','-b':'-b'}
        psym = '0' if p.is_zero() else ('a' if p==AQ else ('-a' if p==-AQ else ('b' if p==BQ else ('-b' if p==-BQ else str(p)))))
        print(f"   pole {psym:>3}: rat  c0={q0[0]}, c1={q0[1]}, c2={q0[2]}")
        print(f"           sqrt2 c0={q1[0]}, c1={q1[1]}, c2={q1[2]}")
