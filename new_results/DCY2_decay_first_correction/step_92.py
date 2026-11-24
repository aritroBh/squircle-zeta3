import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
test = '''
from flint import ctx, arb
ctx.dps = 40
import dcy2_sweep as S
from dcy2_cell import *
for N in [10000, 10100, 10200, 100000, 101000]:
    N = arb(N)
    print(f"N={float(N):.0f}  LogM(-D) = {LogM_of(N,-Delta_of(N)).str(10)}")
Nlo, Nhi = 10000, 10200
A1, R_b = S.deriv_pieces(Nlo, Nhi)
print("A1 =", A1.str(8), " R_b =", R_b.str(8))
print("A0(mid) =", S.A0_of(arb((Nlo+Nhi)/2)).str(10))
'''
open("t_dbg7.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_dbg7.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
