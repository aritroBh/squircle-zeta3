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

H = pickle.load(open('data/towers_D8_K52.pkl','rb'))
H = {k:v for k,v in H.items() if k != '_meta'}
unq4 = kft.unq4
def h_orders(key, k):
    d = key[0]+key[1]
    return {2*d-2-delta: unq4(kd[k]) for delta, kd in H[key].items() if k in kd}

KMAX = 52
# 1) degree-2 kernels vs closed forms
bad20=[]; bad11=[]
for k in range(3, KMAX+1):
    K20, K11, _, _ = kft.kernel_20_21_level(h_orders((3,0),k), h_orders((2,1),k))
    if not kft.pf_equal(K20, kft.R20_closed(k)): bad20.append(k)
    if not kft.pf_equal(K11, kft.R11_closed(k)): bad11.append(k)
print(f"degree-2 kernels vs closed forms k=3..{KMAX}: bad20={bad20}, bad11={bad11}")

# 2) degree-3 kernel law: fit degree-2 polynomials per residue, verify across k=4..KMAX
def res_seqs(nm):
    seqs = {}
    for k in range(4, KMAX+1):
        K30, K21, _, _ = kft.kernel_deg3_level(h_orders((4,0),k), h_orders((3,1),k), h_orders((2,2),k))
        K = K30 if nm=='30' else K21
        assert not [c for c in K[1]+K[2] if not c.is_zero()]
        for p,c in K[0].items():
            seqs.setdefault(p,{})[k]=c
    return seqs

def fit_poly2(ks, ys):
    cur=list(ys); cols=[cur]
    for j in range(1,3):
        cur=[cur[t+1]-cur[t] for t in range(len(cur)-1)]; cols.append(cur)
    k0=ks[0]; d0,d1,d2 = cols[0][0], cols[1][0], cols[2][0]/2
    return (d0 - d1*k0 + d2*Fr(k0*(k0+1)), d1 - d2*Fr(2*k0+1), d2)

def qsi_poly_str(c):
    # c = (c0,c1,c2) Fractions -> string
    c0,c1,c2 = c
    return f"{c2}*k^2 + ({c1})*k + ({c0})"

polesym = {}
for p,name in [(ZERO,'0'),(AQ,'a'),(-AQ,'-a'),(BQ,'b'),(-BQ,'-b')]:
    polesym[p]=name

deg3 = {}
for nm in ['30','21']:
    deg3[nm] = {}
    for p, ck in res_seqs(nm).items():
        ks = sorted(ck); ys = [ck[k] for k in ks]
        # verify degree<=2 across all points
        cur=list(ys); cols=[cur]
        for j in range(1,4):
            cur=[cur[t+1]-cur[t] for t in range(len(cur)-1)]; cols.append(cur)
        le2 = all(c.is_zero() for c in cols[3])
        p0 = fit_poly2(ks, [c.p0 for c in ys])
        p1 = fit_poly2(ks, [c.p1 for c in ys])
        deg3[nm][polesym[p]] = {'rat': [str(x) for x in p0], 'sqrt2': [str(x) for x in p1], 'deg_le_2_all_k': bool(le2)}
print("degree-3 kernel residues all degree<=2 across k=4..52:", all(v['deg_le_2_all_k'] for nm in deg3 for v in deg3[nm].values()))

out = {
 'degree2_kernels_exact_range': f'k=3..{KMAX}',
 'R20_closed': '(1-sqrt2)/(4v) + 1/(4(v-a)) + (2k-7)/(4(v+a)) - (k-2)/(2(v+b))',
 'R11_closed': '(2k-5+sqrt2)/(4v) - (2k-5)/(2(v+a))',
 'degree2_mismatches': {'R20': bad20, 'R11': bad11},
 'degree3_kernel_law_range': f'k=4..{KMAX}',
 'degree3_residue_polynomials': deg3,
}
with open('results/kernel_laws.json','w') as f:
    json.dump(out, f, indent=1)
print("wrote results/kernel_laws.json")
with open('progress.log','a') as f:
    f.write(f"step_20 Task C with K=52: deg-2 kernels==closed forms k=3..{KMAX} (bad20={bad20},bad11={bad11}); deg-3 kernel law k=4..{KMAX}\n")
