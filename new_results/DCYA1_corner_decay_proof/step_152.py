import sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
from arb_ei import E_I_cell, E_I_sc, iv, arb, alpha
Ncell=iv("24468","26057")
EI=E_I_cell(Ncell)
n2=arb(Ncell).upper()
sc=n2**arb("-1.5")*(-(2*alpha*n2).sqrt()).exp()
print("E_I:",EI)
print("sc:",sc)
print("E_I/sc:",arb(EI/sc))
# point value at N=25000 via thin cell
print("thin [25000,25010]:", E_I_sc(iv("25000","25010")))
