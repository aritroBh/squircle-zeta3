import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
test = '''
from flint import ctx, arb
ctx.dps = 40
from dcy2_cell import *
for N in [83966.56, 100000.0]:
    print(f"N={N}")
    print("  y =", y_of(N).str(10), " Delta =", Delta_of(N).str(10), " kap =", kap_of(N).str(8))
    print("  alpha =", alpha_of(N).str(10))
    print("  L3 =", L_of(N,3).str(10), " L4*D =", (L_of(N,4)*Delta_of(N)).str(10), " L5*D^2 =", (L_of(N,5)*Delta_of(N)**2).str(10))
    print("  al*D/(c/2) =", (alpha_of(N)*Delta_of(N)/(c2min_of(N)/2)).str(10))
    print("  B3*D/(c/2) =", (B3_of(N)*Delta_of(N)/(c2min_of(N)/2)).str(10))
'''
open("t_dbg9.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_dbg9.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
