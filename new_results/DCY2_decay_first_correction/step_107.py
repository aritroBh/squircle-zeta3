import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
test = '''
from flint import ctx, arb
ctx.dps = 40
from dcy2_cell import *
from dcy2_sweep import TL_cell, TR_cell, lamL_cell
import dcy2_cert as CC
for Nf in [400, 10000, 100000]:
    N = arb(Nf)
    q = CC.quantities(N)
    print(f"N={Nf}: budget Tails = {q['Tails'].str(6)}")
    hi = Nf*1.001
    tl = TL_cell(Nf, hi)
    tr = TR_cell(Nf, hi)
    lhs_tl = (2/PI)*y_of(N)*S_sup(N)*tl*N**2
    lhs_tr = (2/PI)*y_of(N)*S_sup(N)*tr*N**2
    print(f"   sweep TL_cell = {tl.str(6)}, TR_cell = {tr.str(6)}")
    print(f"   LHS contrib TL = {lhs_tl.str(8)}, TR = {lhs_tr.str(8)}")
    print(f"   lamL = {lamL_of(N).str(8)}  LogM(-D) = {LogM_of(N,-Delta_of(N)).str(8)}")
'''
open("t_dbg13.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_dbg13.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
