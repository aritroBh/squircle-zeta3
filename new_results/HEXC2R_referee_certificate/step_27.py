
import subprocess
code = '''
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../HEXC1_hexagon_certificate")
import cmzv6, reflib as R
lower_m, lower_t = {}, {}
for n in (1,2,3):
    rt, ct = cmzv6.make_relations(n, lower_t)
    rm, cm = R.make_relations(n, lower_m)
    lower_t[n] = rt; lower_m[n] = rm
print(len(rt), len(rm), len(set(rt)-set(rm)), len(set(rm)-set(rt)))
'''
for i in range(3):
    r = subprocess.run(["<project>/.venv/bin/python", "-c", code], capture_output=True, text=True)
    print("run", i, ":", r.stdout.strip(), r.stderr[-300:] if r.returncode else "")
