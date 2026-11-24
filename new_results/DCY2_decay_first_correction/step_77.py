import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
test = '''
from flint import ctx, arb
ctx.dps = 40
import dcy2_cert as C
N = arb(460, 60)
q = C.quantities(N)
for k, v in q.items():
    try:
        print(k, "=", v.str(10) if hasattr(v, "str") else v)
    except Exception as e:
        print(k, "ERR", e)
'''
open("t_dbg2.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_dbg2.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
