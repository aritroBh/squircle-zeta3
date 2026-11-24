import os, sys, subprocess, datetime
def log(msg):
    with open('progress.log','a') as fh:
        fh.write(f"[{datetime.datetime.now().isoformat()}] {msg}\n")

# check P.num exists
src = open('polylog_exact.py').read()
print('def num' in src)
i = src.find('def num')
print(src[i:i+700])
log("omega4.json saved: exact (pole,residue) data for all k=4 pieces; reconstruction asserted exact in sympy.")
