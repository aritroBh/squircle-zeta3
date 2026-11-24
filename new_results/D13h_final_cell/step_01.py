import shutil, subprocess, os, sys

SRC = 'research/08_theory/D13g_all_N/step_38.py'
WORK = os.getcwd()
print("CWD:", WORK)
dst = os.path.join(WORK, 'step_38.py')
shutil.copy(SRC, dst)
print("copied to", dst)

# Run its __main__
r = subprocess.run([sys.executable, dst], capture_output=True, text=True, timeout=120)
print("STDOUT:")
print(r.stdout)
print("STDERR:")
print(r.stderr[-2000:] if r.stderr else "(none)")
