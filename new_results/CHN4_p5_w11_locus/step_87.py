import subprocess, time
def log_progress(line):
    with open('progress.log','a') as f: f.write(time.strftime('%H:%M:%S ')+line+'\n')
log_progress("det11 rank re-verified = 693 mod 2 primes. Divisibility proof airtight. Launching control: zeta3*zeta(3,5) coordinate (W11 idx7) interpolation + divisibility by a2a3a4(Sum a) and C4.")
proc=subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','run_ctrl_z3z35.py'],
                      stdout=open('run_ctrl_z3z35.out','w'), stderr=subprocess.STDOUT)
open('run_ctrl_z3z35.pid','w').write(str(proc.pid))
print("launched control pid", proc.pid)
