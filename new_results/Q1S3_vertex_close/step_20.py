import shutil, os
src = "<project>/research/08_theory/Q1S_vertex_symbol/"
for f in ["q1s_exact.py", "q1s_symbol.py", "sunit_li3.pkl"]:
    shutil.copy(src+f, ".")
print(os.listdir('.'))
# quick self-test of the exact engine
import subprocess
print(subprocess.run(["python3","q1s_exact.py"],capture_output=True,text=True).stdout)
