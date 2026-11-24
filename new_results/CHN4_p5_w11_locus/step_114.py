import time
def log_progress(line):
    with open('progress.log','a') as f: f.write(time.strftime('%H:%M:%S ')+line+'\n')
log_progress("Final confirmations: P7 irreducible (factor_list), reversal-symmetric (monomial check 0/209 + subs simultaneous=True), P7.txt==-10*Q11, f11=-(1/10)a2a3a4(Sum a)P7 verified vs engine at 3 random pts. NOTE: sympy subs({x:y,y:x}) is SEQUENTIAL in this build (trap); use simultaneous=True or monomial-level checks for reversal symmetry. C4,Q7' also reversal-symmetric. summary.md finalized.")
print("done")
import os
print("work dir files:")
for f in sorted(os.listdir('.')):
    if f.endswith(('.txt','.pkl','.md','.log')) and not f.startswith('step'):
        print("  ", f, os.path.getsize(f))
