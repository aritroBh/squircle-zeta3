import subprocess, sys
r = subprocess.run([sys.executable, "step_04.py"], capture_output=True, text=True, cwd="research/08_theory/OCT7_generic_octagon")
print(r.stdout)
print("STDERR:", r.stderr[-500:] if r.stderr else "none")
