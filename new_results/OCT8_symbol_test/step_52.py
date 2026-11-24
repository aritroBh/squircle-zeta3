import mpmath as mp, pickle, sympy as sp, re, sys
from fractions import Fraction as Fr
sys.path.insert(0,".")
mp.mp.dps = 300
def M(q):
    q = Fr(q); return mp.mpf(q.numerator)/q.denominator
NS = pickle.load(open("nullspace.pkl","rb"))["NS"]
dB = pickle.load(open("symbol_identities_B.pkl","rb"))
G_names, D_names = dB["G_names"], dB["D_names"]
logv = {p: mp.log(p) for p in (2,3,5)}
Li3v = [mp.polylog(3, M(Fr(gn[4:-1]))) for gn in G_names]
Dv = []
for dn in D_names:
    m = re.match(r"-Li2\((.*)\)\*log(\d+)", dn)
    w, p = Fr(m.group(1)), int(m.group(2))
    Dv.append(-mp.polylog(2, M(w))*logv[p])
vals = []
for nv in NS:
    tot = mp.mpc(0)
    for i in range(len(G_names)):
        if nv[i] != '0': tot += M(Fr(nv[i]))*Li3v[i]
    for j in range(len(D_names)):
        if nv[len(G_names)+j] != '0': tot += M(Fr(nv[len(G_names)+j]))*Dv[j]
    vals.append(tot)
# real and imaginary parts
z3, z2 = mp.zeta(3), mp.zeta(2)
# build independent real value basis greedily
basis_vals = []   # list of (value, description)
def try_add(v, desc):
    for b,_ in basis_vals:
        # check independence crudely via pslq
        pass
    basis_vals.append((v, desc))
rev = [mp.re(v) for v in vals]
imv = [mp.im(v) for v in vals]
# find distinct nonzero real values
import itertools
nonzero = [v for v in rev if abs(v) > mp.mpf('1e-250')]
print("#null vecs:", len(vals), " #nonzero real values:", len(nonzero))
# cluster: PSLQ each against {z3, z2log p}
inb = [z3, z2*logv[2], z2*logv[3], z2*logv[5]]
out_vals = []
for v in nonzero:
    rel = mp.pslq([v]+inb, tol=mp.mpf('1e-260'), maxcoeff=10**60, maxsteps=1000000)
    if rel is None:
        out_vals.append(v)
print("values outside {z3,z2log}:", len(out_vals))
for v in out_vals[:8]:
    print("   ", mp.nstr(v, 18))
# imaginary parts span
nz_im = [v for v in imv if abs(v) > mp.mpf('1e-250')]
print("#nonzero imaginary:", len(nz_im))
for v in nz_im[:5]:
    print("   im:", mp.nstr(v, 18))
