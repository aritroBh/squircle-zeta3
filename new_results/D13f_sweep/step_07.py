
import subprocess, os
os.chdir("research/08_theory/D13f_sweep")
PY = "<home>)/.venv/bin/python"
# sanity
r = subprocess.run([PY, "-c", "import flint; from flint import arb; print('flint', flint.__version__ if hasattr(flint,'__version__') else 'ok'); print(arb.pi())"],
                   capture_output=True, text=True)
print("STDOUT:", r.stdout)
print("STDERR:", r.stderr[:2000])
