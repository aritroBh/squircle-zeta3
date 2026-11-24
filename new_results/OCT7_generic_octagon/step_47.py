import subprocess, sys
for s in ["step_04.py","step_05.py"]:
    r = subprocess.run([sys.executable, s], capture_output=True, text=True, cwd="research/08_theory/OCT7_generic_octagon")
    print("="*20, s, "="*20)
    print(r.stdout[-1800:])
    if r.stderr: print("STDERR:", r.stderr[-300:])
