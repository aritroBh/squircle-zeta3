import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_05f.py"],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=1500)
print(r.stdout); print("ERR:", (r.stderr or "")[-1000:])
