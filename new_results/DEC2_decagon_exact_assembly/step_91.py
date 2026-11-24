import subprocess, time
t0=time.time()
r = subprocess.run(['<project>/.venv/bin/python',
                    '<project>/research/08_theory/PRV1_depth2_proofs/verify_full.py'],
                   capture_output=True, text=True, timeout=280)
print("runtime %.0fs"%(time.time()-t0))
print(r.stdout[-1500:])
print("STDERR tail:", r.stderr[-300:] if r.stderr else "")
