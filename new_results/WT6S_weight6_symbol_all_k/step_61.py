
import subprocess, time, os
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
PY = '<project>/.venv/bin/python'
for st in ['step_04.py', 'step_05.py', 'step_06.py']:
    t0 = time.time()
    r = subprocess.run(['nice', '-n', '10', PY, '-B', st], capture_output=True, text=True, cwd=WD)
    print(f"{st}: {time.time()-t0:.1f}s  rc={r.returncode}")
