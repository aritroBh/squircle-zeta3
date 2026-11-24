import sys, subprocess
cor = "<project>/research/08_theory/COR4_M3alpha_classical/"
r = subprocess.run([sys.executable, cor+"step12_specializations.py"], capture_output=True, text=True, timeout=300)
print(r.stdout)
print(r.stderr[-800:] if r.returncode else "exit 0")
