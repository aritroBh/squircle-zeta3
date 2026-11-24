
import os, shutil

BASE = 'research/08_theory'
SRC = os.path.join(BASE, 'D13f_sweep')
WORK = os.path.join(os.getcwd())
print("CWD:", WORK)
print("D13f_sweep contents:", os.listdir(SRC))
print("D13d contents:", os.listdir(os.path.join(BASE, 'D13d_decay_theorem')))

# Copy the files I'm allowed to read
for f in ['step_37.py', 'step_34.py', 'result.md']:
    p = os.path.join(SRC, f)
    if os.path.exists(p):
        shutil.copy(p, os.path.join(WORK, f))
        print("copied", f, os.path.getsize(p))
