import os
os.makedirs("research/08_theory/COR6X_independent_corner_symbols", exist_ok=True)
os.chdir("<project>/research/08_theory/COR6X_independent_corner_symbols")
open("progress.log","a").write("start COR6X referee re-implementation\n")
open("result.md","w").write("# COR6X independent referee report (draft)\n\nWork in progress.\n")
print(open("step_01.py").read() if os.path.exists("step_01.py") else "no step_01")
