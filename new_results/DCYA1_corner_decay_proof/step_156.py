import sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
from flint import arb, acb, ctx
ctx.dps=60
from arb_ei import E_I_cell, E_I_sc, iv, alpha
Ncell=iv("20000","20142")
EI=E_I_cell(Ncell)
print("E_I_cell:",EI)
n2=arb(Ncell).upper(); sc=n2**arb("-1.5")*(-(2*alpha*n2).sqrt()).exp()
print("E_I/sc:",arb(EI/sc))
print("E_I_sc:",E_I_sc(Ncell))
# also at [20000,20050]
EI2=E_I_cell(iv("20000","20050"))
sc2=arb("20050")**arb("-1.5")*(-(2*alpha*arb("20050")).sqrt()).exp()
print("E_I/sc [20000,20050]:",arb(EI2/sc2))
