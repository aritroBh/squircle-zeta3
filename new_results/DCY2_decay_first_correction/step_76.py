import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
test = '''
from flint import ctx, arb
ctx.dps = 40
import dcy2_cert as C
def cell(lo, hi):
    N = arb((lo+hi)/2, (hi-lo)/2)
    r, q = C.ratio(N)
    return r
for lo, hi in [(400, 520), (520, 676), (1000, 1300), (2000, 2600), (10000, 13000), (10**5, 1.3e5), (10**6, 1.3e6)]:
    r = cell(lo, hi)
    print(f"[{lo:9d}, {hi:9d}]  ratio upper = {r.upper().str(12)}")
'''
open("t_cells.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_cells.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
