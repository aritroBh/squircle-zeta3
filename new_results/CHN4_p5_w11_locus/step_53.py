import os, time
# check background proof
print("=== proof_div.log ===")
print(open('proof_div.log').read() if os.path.exists('proof_div.log') else "not yet")
def log_progress(line):
    with open('progress.log','a') as f: f.write(time.strftime('%H:%M:%S ')+line+'\n')
log_progress("Stage3 sign survey: built robust high-precision V10 solver (exact a4 isolation via G.intervals + mp bisection + common a5 root). CONFIRMED at 60 digits: Q11 CHANGES SIGN on V10+ : Q11=-0.0495 at (1/100,1/100), Q11=-140429 at (1,1), Q11=+0.222 at (1/100,3). Also 120x120 numpy scan over [0.01,100]^2: ~12000 pos + ~11700 neg interior pts. So Gamma=V10+ cap {f11=0} is a REAL curve -> go to Step4 (W12). NOTE: a5=-B/A gives SPURIOUS lifts at degenerate a4 (e.g. a4=1 double root gives fake a5=0.344); must use common-root a5.")
print("progress updated")
