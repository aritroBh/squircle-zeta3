import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
test = '''
from flint import ctx, arb
ctx.dps = 40
import dcy2_sweep as S
from dcy2_cell import *
for lo, hi in [(400, 408), (2000, 2040)]:
    N = arb((lo+hi)/2, (hi-lo)/2)
    em = S.E_mom_iv(N)
    et = S.E_tailmom_iv(N)
    print(f"[{lo},{hi}]")
    print("  E_mom  mid/rad:", em.mid(), em.rad())
    print("  E_tail mid/rad:", et.mid(), et.rad())
    v = S.LHS_cell_A(lo, hi)
    print("  LHS mid/rad:", v.mid(), v.rad())
'''
open("t_midrad.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_midrad.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
