import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
test = '''
from flint import ctx, arb
ctx.dps = 40
from dcy2_cell import *
from dcy2_sweep import *
Nlo, Nhi = 100000, 102000
br = window_cell(Nlo, Nhi)
print("window bracket:", br.str(10))
tp = TR_cell(Nlo, Nhi) + TL_cell(Nlo, Nhi)
print("tailpieces:", tp.str(10))
print("E_mom alone:", E_mom_of(arb(Nlo)).str(10), " E_tailmom:", E_tailmom_of(arb(Nlo)).str(10))
# pieces of window_cell
y_lo = y_of(Nlo)
al = alpha_of(Nlo); B3 = B3_of(Nlo)
print("alpha:", al.str(8), " B3:", B3.str(8))
c_lo, c_hi = c2min_of(Nlo), c2min_of(Nhi)
D_lo, D_hi = Delta_of(Nlo), Delta_of(Nhi)
print("al*D_lo/(c_lo/2) =", (al*D_lo/(c_lo/2)).str(8))
print("B3*D_lo/(c_lo/2) =", (B3*D_lo/(c_lo/2)).str(8))
print("Mtilde_of(Nhi, 12, c_hi/2):", Mtilde_of(Nhi, 12, c_hi/2).str(8))
print("al^4/24 * GsplitB(12,al):", (al**4/24 * (Mtilde_of(Nhi,12,c_hi/2) if al*D_lo <= c_lo/2 else -1)).str(8))
print("M4p:", M4p_of(Nlo).str(8))
print("GsplitB(4,B3) =", (Mtilde_of(Nhi,4,c_hi/2) if B3*D_lo <= c_lo/2 else -1).str(8))
print("H_up:", H_up_of(Nlo).str(8))
print("S_sup:", S_sup(Nlo).str(8))
print("T0:", T0_cell(Nlo,Nhi).str(8))
'''
open("t_dbg4.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_dbg4.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
