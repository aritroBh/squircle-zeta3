import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
import mpmath as mp
from flint import fmpz_mat
from oct5d3_basis import atom_specs, C_SPECS
from oct5d3_atoms import INVISIBLE5, prod, k32_value

sf = pickle.load(open("symfit_m4f.pkl","rb"))
x0, labels = sf['x0'], sf['vis_labels']
Kred = pickle.load(open("kernel_red.pkl","rb"))['Kred']
mixed = pickle.load(open("mixed61_vals.pkl","rb"))
spec_map = {lab: spec for lab, kind, spec in atom_specs()}
vis = [(lab, spec_map[lab]) for lab in labels]
def q2f(z):
    return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*mp.sqrt(2)

D = 1500
mp.mp.dps = D
A = {k: mp.mpf(v) for k, v in pickle.load(open("atoms_1650.pkl","rb")).items()}
V = [prod(A, spec) for lab, spec in vis]
r2 = mp.sqrt(2)
kvals = []
for row in Kred:
    v = mp.mpf(0)
    for j in range(len(labels)):
        p, q = row[2*j], row[2*j+1]
        if p or q: v += (mp.mpf(p) + mp.mpf(q)*r2)*V[j]
    kvals.append(v)
atom_names = [nm for nm, sp in INVISIBLE5] + ["K32"] + [f"MIX{i+1}" for i in range(7)] + [f"K{i+1}" for i in range(22)]
atom_vals = ([prod(A, sp) for nm, sp in INVISIBLE5] + [k32_value()]
             + [mp.mpf(val) for val, desc in mixed] + kvals)
print("51 atoms evaluated")

# build the kernel lattice over Q(sqrt2): rows [e_i | C*w_i], w alternates (v, r2*v)
t0 = time.time()
wvals = []
for v in atom_vals: wvals += [v, r2*v]
n = len(wvals)
C = mp.mpf(10)**int(D*0.92)
rows = []
for i in range(n):
    row = [0]*(n+1); row[i] = 1; row[n] = int(mp.nint(C*wvals[i])); rows.append(row)
M = fmpz_mat(rows).lll()
print("LLL done %.0fs" % (time.time()-t0))
# short rows with last coord 0 -> relations
rels = []
for i in range(n):
    last = int(M[i, n])
    v = [int(M[i,j]) for j in range(n)]
    h = max(abs(x) for x in v)
    res = abs(sum(mp.mpf(v[j])*wvals[j] for j in range(n)))
    if res < mp.mpf(10)**(-D+150) and h < 10**30:
        rels.append((h, res, v))
print("candidate atom-relations:", len(rels))
for h, res, v in rels[:12]:
    terms = [(atom_names[j//2], (v[j] if j%2==0 else 0, v[j] if j%2==1 else 0)) for j in range(n) if v[j]]
    terms = [(nm, (v[2*k], v[2*k+1])) for k, (nm, _) in enumerate([(atom_names[j//2], None) for j in range(n) if v[j]])] if False else None
    print("  h=%s res=%s nterms=%d" % (mp.nstr(mp.mpf(h),2), mp.nstr(res,2), sum(1 for x in v if x)))
pickle.dump({'rels': [(h, [int(x) for x in v]) for h, res, v in rels], 'names': atom_names}, open("atom_kernel.pkl","wb"))
print("saved atom_kernel.pkl")
