import os, datetime
print(os.getcwd())
print(os.listdir('.'))
BASE = '<project>'
WORK = os.getcwd()

def log(msg):
    with open('progress.log','a') as f:
        f.write(f"[{datetime.datetime.now().isoformat()}] {msg}\n")

log("Starting M4a tasks in existing work dir.")

with open(os.path.join(BASE,'paper/squircle_note/main.tex')) as f:
    tex = f.read()
print('tex length', len(tex))
import re
for key in ['prop:split','rem:k4','sec:M3','eq:W4','W4']:
    idxs = [m.start() for m in re.finditer(re.escape(key), tex)]
    print(key, idxs[:12])
