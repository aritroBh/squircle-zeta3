import os, datetime
BASE = '<project>'
WORK = os.path.join(BASE, 'work_m4a')
os.makedirs(WORK, exist_ok=True)
os.chdir(WORK)

# progress.log
def log(msg):
    with open('progress.log','a') as f:
        f.write(f"[{datetime.datetime.now().isoformat()}] {msg}\n")

log("Work directory created. Starting M4a tasks.")

# Read main.tex
with open(os.path.join(BASE,'paper/squircle_note/main.tex')) as f:
    tex = f.read()
print(len(tex))
# Find prop:split and W4 and rem:k4
import re
# print sections around these labels
for key in ['prop:split','rem:k4','W4','sec:M3']:
    idxs = [m.start() for m in re.finditer(re.escape(key), tex)]
    print(key, idxs[:10])
