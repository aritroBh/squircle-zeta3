
import sys, os, subprocess
print(sys.executable)
print(sys.version)
import sympy, mpmath
print("sympy", sympy.__version__, "mpmath", mpmath.__version__)
base = "<project>"
# check venv python
r = subprocess.run([base + "/.venv/bin/python", "-c", "import sys, sympy, mpmath; print(sys.executable); print(sympy.__version__, mpmath.__version__)"], capture_output=True, text=True)
print(r.stdout, r.stderr)
# logs of TOW4
tow4 = base + "/research/08_theory/TOW4_tower_depth8"
for f in ["verify_tow4.log", "forced78.log", "check_full78.log", "adj78_all.log"]:
    print("="*20, f)
    print(open(tow4 + "/" + f).read())
