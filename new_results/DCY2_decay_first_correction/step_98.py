import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
test = '''
from flint import ctx, arb
ctx.dps = 40
import dcy2_sweep as S
from dcy2_cell import *
lo, hi = 83966.56115454846, 85645.89237763942
br = S.window_cell(lo, hi)
tp = S.TR_cell(lo,hi) + S.TL_cell(lo,hi)
print("window bracket:", br.str(10))
print("tailpieces:", tp.str(10))
y = y_of(lo)
print("prefactor:", ((2/PI)*y*S_sup(lo)*H_up_of(lo)).str(8))
print("H_up*bracket:", (H_up_of(lo)*br).str(8))
print("pointwise window at 84000:", (E_mom_of(arb(84000)) + E_tailmom_of(arb(84000))).str(10))
# E_mom pieces at Nlo
al = alpha_of(lo); B3 = B3_of(lo)
c_lo, c_hi = c2min_of(lo), c2min_of(hi)
D_lo = Delta_of(lo)
print("al*D_lo/(c_lo/2) =", (al*D_lo/(c_lo/2)).str(8))
print("B3*D_lo/(c_lo/2) =", (B3*D_lo/(c_lo/2)).str(8))
print("Mtilde(Nhi,12,c_hi/2) =", Mtilde_of(hi,12,c_hi/2).str(8))
print("Mtilde(Nhi,4,c_hi/2) =", Mtilde_of(hi,4,c_hi/2).str(8))
print("M4p(lo) =", M4p_of(lo).str(8))
'''
open("t_dbg8.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_dbg8.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
