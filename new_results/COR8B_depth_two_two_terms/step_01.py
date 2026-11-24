
import os, time
os.makedirs('steps', exist_ok=True)
open('progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_00 start: plan + workspace setup\n")
# quick sanity of available libs
import sympy, mpmath, numpy, scipy
print("sympy", sympy.__version__, "| mpmath", mpmath.__version__)
try:
    import flint
    print("flint", flint.__version__)
except Exception as e:
    print("flint import failed:", e)
