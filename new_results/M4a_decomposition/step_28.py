import subprocess, sys, datetime
def log(msg):
    with open('progress.log','a') as fh:
        fh.write(f"[{datetime.datetime.now().isoformat()}] {msg}\n")

# Cross-check: re-run ibp_reduce.py 4 in my folder
r = subprocess.run([sys.executable, 'ibp_reduce.py', '4'], capture_output=True, text=True, timeout=1200)
open('ibp_reduce_k4.out','w').write(r.stdout)
print(r.stdout)
print('STDERR:', r.stderr[:500])
log("Ran ibp_reduce.py 4 -> ibp_reduce_k4.out (cross-check of paper output).")
