import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
test = '''
from flint import ctx
ctx.dps = 40
import dcy2_cert as C
for N in [2000, 10000, 10**6]:
    q = C.quantities(N)
    print(f"N={N}")
    for k in ["E_mom","E_tailmom","Tails","eRe_phi","H_up","V3","VA","Delta","sig","kap"]:
        print("   ", k, q[k].str(8) if hasattr(q[k],"str") else q[k])
    br = (2/C.PI)*q["eRe_phi"]*q["H_up"]
    print("    prefactor*(E_mom+E_tailmom) =", (br*(q["E_mom"]+q["E_tailmom"])).str(8))
    print("    budget =", C.budget(N).str(8))
'''
open("t_debug.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_debug.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
