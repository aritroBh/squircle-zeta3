import subprocess, time
t0 = time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python',
                    '<project>/research/08_theory/DEC2_decagon_exact_assembly/verify_dec2.py'],
                   capture_output=True, text=True, timeout=295)
print("runtime %.0fs, returncode %d" % (time.time()-t0, r.returncode))
print(r.stdout)
if r.stderr: print("STDERR tail:", r.stderr[-2000:])
