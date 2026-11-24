import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
import mpmath as mp
from oct5d3_atoms import W1_BASIS, prod
from q2_lll import fit
from flint import fmpz_mat

mp.mp.dps = 800
Afull = {k: mp.mpf(v) for k, v in pickle.load(open("atoms_1650.pkl","rb")).items()}
Wt = pickle.load(open("wparts_true_1650.pkl","rb"))
T = mp.mpf(Wt['Wt1'])
atoms = [(nm, prod(Afull, sp)) for nm, sp in W1_BASIS]
print("T =", mp.nstr(T, 30))
print("atoms:", [(nm, mp.nstr(v, 20)) for nm, v in atoms])
# residual of the known relation with THESE values
r2 = mp.sqrt(2)
res = T - ((4800-3360*r2)*atoms[0][1] + (2880-1920*r2)*atoms[1][1])
print("known-relation residual with these values:", mp.nstr(abs(res), 3))
# build the lattice the way fit() does and print rows
vals = [T]
for nm, v in atoms: vals += [v, r2*v]
n = len(vals)
C = mp.mpf(10)**736
rows = []
for i in range(n):
    row = [0]*(n+1); row[i] = 1; row[n] = int(mp.nint(C*vals[i])); rows.append(row)
M = fmpz_mat(rows).lll()
print("LLL row0:", [int(M[0,j]) for j in range(n)])
