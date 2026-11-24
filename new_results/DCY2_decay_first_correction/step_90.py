import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
test = '''
from flint import ctx, arb
ctx.dps = 40
import dcy2_sweep as S
from dcy2_cell import *
Nlo, Nhi = 100000, 102000
print("TR_cell:", S.TR_cell(Nlo, Nhi).str(10))
print("lamL_cell:", S.lamL_cell(Nlo, Nhi).str(10))
Nmid = (Nlo+Nhi)/2; hw = (Nhi-Nlo)/2
V = S.dLam_dN(arb(Nmid, hw))
print("V =", V.str(10))
LM = S.LogM_mD_mv(Nmid, V, hw)
print("LM =", LM.str(10))
CH_L = (PI/2).sqrt()*2**(arb(3)/4)
print("TL pieces: exp(LM) =", LM.exp().str(8), " CH_L/lamL =", (CH_L/S.lamL_cell(Nlo,Nhi)).str(8))
print("TL_cell:", S.TL_cell(Nlo, Nhi).str(10))
'''
open("t_dbg6.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_dbg6.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
