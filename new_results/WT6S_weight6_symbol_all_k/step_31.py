
import os, sys, subprocess
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
os.chdir(WD)
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice', '-n', '10', PY, '-B', 'wt6s_core.py'], capture_output=True, text=True, cwd=WD)
print(r.stdout)
print("STDERR:", r.stderr[-2000:])
