
import os, json, subprocess
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
os.chdir(wd)
# Check existing step_01.py
print(open("step_01.py").read()[:1500])
print("=== progress.log exists:", os.path.exists("progress.log"))
PY = "<project>/.venv/bin/python"
print(subprocess.run([PY,"-c","import mpmath, sympy, flint; print(mpmath.__version__, sympy.__version__, flint.__version__)"],capture_output=True,text=True).stdout)
