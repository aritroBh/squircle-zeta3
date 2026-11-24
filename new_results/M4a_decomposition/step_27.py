import os, shutil, datetime
BASE = '<project>'
SCR = os.path.join(BASE,'paper/squircle_note/scripts')
shutil.copy(os.path.join(SCR,'dilog_exact.py'), 'dilog_exact.py')

def log(msg):
    with open('progress.log','a') as fh:
        fh.write(f"[{datetime.datetime.now().isoformat()}] {msg}\n")
log("Copied dilog_exact.py too (polylog_exact dependency).")
print(open('dilog_exact.py').read()[:800])
