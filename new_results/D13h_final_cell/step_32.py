import subprocess, sys
# Final end-to-end sanity: fresh runs of the three key certificates
for f, key in [('step_10.py','FINAL CELL'), ('step_15.py','sup C'), ('step_16.py','match')]:
    r = subprocess.run([sys.executable, f], capture_output=True, text=True, timeout=600)
    lines = [l for l in r.stdout.splitlines() if key in l or 'FINAL CELL' in l]
    print(f"--- {f} ---")
    for l in lines[:4]: print("  ", l)
    if r.stderr: print("  STDERR:", r.stderr[-200:])
print("\nprogress.log tail:")
print(open('progress.log').read()[-500:])
