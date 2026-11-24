import os, shutil, datetime
BASE = '<project>'
SCR = os.path.join(BASE,'paper/squircle_note/scripts')
WORK = os.getcwd()
for f in ['mk_integrands.py','ibp_reduce.py','qsqrt2.py','polylog_exact.py','m3_pslq.py']:
    shutil.copy(os.path.join(SCR,f), os.path.join(WORK,f))
print(sorted(os.listdir(WORK)))

def log(msg):
    with open('progress.log','a') as fh:
        fh.write(f"[{datetime.datetime.now().isoformat()}] {msg}\n")
log("Copied paper scripts into work dir; read all reference material (main.tex, w4_symbol, cmzv, verify_D, check_conversion, A3b results).")

# quick peek at polylog_exact P class
src = open('polylog_exact.py').read()
print(src[:3000])
