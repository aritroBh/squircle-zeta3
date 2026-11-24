import subprocess, os, json, pickle
VENV = "<project>/.venv/bin/python"
r = subprocess.run([VENV, "-c", "import sys, sympy, mpmath, flint; print(sys.version); print('sympy', sympy.__version__, 'mpmath', mpmath.__version__, 'flint ok')"], capture_output=True, text=True)
print(r.stdout, r.stderr)
base = "<project>/research/08_theory/G2_M4alpha/"
print(open(base+"b4_pi4_part.json").read())
