import subprocess, time, os
T0 = time.monotonic()
base = '<project>/research/08_theory'
py = '<project>/.venv/bin/python'
r = subprocess.run(['nice', py, 'verify_cor13c.py'], cwd=os.path.join(base,'COR13C_third_corner_certificate'),
                   capture_output=True, text=True, timeout=1200)
print(r.stdout[-3500:])
print("rc:", r.returncode, f"({time.monotonic()-T0:.0f}s)")
with open('progress.log','a') as f:
    f.write("referee step(their verifier): fresh verify_cor13c.py run: rc=%d (%.0fs)\n" % (r.returncode, time.monotonic()-T0))
