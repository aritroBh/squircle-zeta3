import subprocess, time
# record the W12 gamma results
def log_progress(line):
    with open('progress.log','a') as f: f.write(time.strftime('%H:%M:%S ')+line+'\n')
log_progress("Stage4 W12 batch DONE (6 Gamma pts, a1=0.0005..0.2): at EVERY pt f11 & z3z35 vanish (~1e-49..1e-72) AND all W12 irreducibles NONZERO (z2^2z35,z2z37,CLASS,z2118). max|W12 irred| ranges 35.5 (a1=0.0005) to 70757 (a1=0.1). => NO Gamma pt is irreducible-free through W12. ANSWER: NO positive point is irreducible-free through W12 (numerical).")
# launch divisibility proof (heavy, background)
proc=subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','run_proof_div.py'],
                      stdout=open('run_proof_div.out','w'), stderr=subprocess.STDOUT)
open('run_proof_div.pid','w').write(str(proc.pid))
print("launched divisibility proof pid", proc.pid)
log_progress("Launched 693-point rigorous divisibility proof (D=f11-L*Q11 on rank-693 reversal-inv deg-11 det set) in background.")
