import subprocess, sys, os
base = '<project>/research/08_theory'
# find w5pipe
r = subprocess.run(['find', f'{base}/ZQ5b_zagier_w5_larger_S', '-name', 'w5pipe*'], capture_output=True, text=True)
print(r.stdout, r.stderr)
# Check python version and packages
print("python:", sys.executable, sys.version)
try:
    import mpmath; print("mpmath", mpmath.__version__)
except Exception as e: print("mpmath ERR", e)
try:
    import flint; print("flint", flint.__version__)
except Exception as e: print("flint ERR", e)
try:
    import sympy; print("sympy", sympy.__version__)
except Exception as e: print("sympy ERR", e)
