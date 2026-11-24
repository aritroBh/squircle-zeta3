src = open("mixed_engine_gen.py").read()
old = 'log_23 = l3 - l2     # log(2/3) = log((k1-k2)/k1), as symbols'
new = 'log_23 = l2 - l3     # log(2/3) = log((k1-k2)/k1), as symbols'
assert old in src; src = src.replace(old, new)
open("mixed_engine_gen.py","w").write(src)
import subprocess, sys
for s in ["step_04.py","step_05.py"]:
    r = subprocess.run([sys.executable, s], capture_output=True, text=True, cwd="research/08_theory/OCT7_generic_octagon")
    out = r.stdout
    print("="*20, s)
    for line in out.splitlines():
        if any(k in line for k in ["worst","residual","diff:","E_arc (exact","independent quad","Koelbig"]):
            print(line)
