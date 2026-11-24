import subprocess, os, shutil
cwd='<project>/research/08_theory/DCYA1_corner_decay_proof/'
os.chdir(cwd)
# Initialize progress.log and result.md
with open('progress.log','w') as f:
    f.write("# DCYA1 progress log\n")
with open('result.md','w') as f:
    f.write("# DCYA1 — Decay law of M_k(alpha) at every corner angle (PROOF)\n\n(work in progress)\n")
# Copy corner_series.py reference for import (read-only source; copy logic into workspace)
shutil.copy('<project>/research/10_frontier/G1_zeta3_any_angle/corner_series.py', cwd+'corner_series.py')
print("Workspace set up. Files:", os.listdir(cwd))
print("corner_series copied.")
PY='<project>/.venv/bin/python'
print("python:", PY)
print(subprocess.run([PY,'-c','import mpmath, flint; print("mpmath",mpmath.__version__); print("flint ok")'], capture_output=True, text=True).stdout)
