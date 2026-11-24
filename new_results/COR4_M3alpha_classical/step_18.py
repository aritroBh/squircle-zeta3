import sys, subprocess
print(sys.executable)
cor = "<project>/research/08_theory/COR4_M3alpha_classical/"
r = subprocess.run([sys.executable, cor+"step10_symbol_checks.py"], capture_output=True, text=True)
print(r.stdout)
print(r.stderr[-1500:] if r.returncode else "exit 0")
