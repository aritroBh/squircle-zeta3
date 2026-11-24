import subprocess, os, time
wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
code = r'''
import flint
help(flint.acb.integral)
'''
r = subprocess.run(["<project>/.venv/bin/python","-c", code], capture_output=True, text=True, cwd=wd)
print(r.stdout[:2500]); print(r.stderr[:300])
