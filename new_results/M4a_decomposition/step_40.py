import subprocess, sys
r = subprocess.run([sys.executable, 'm3_tails_exact.py', '3'], capture_output=True, text=True, timeout=2400)
tail = r.stdout[-1800:]
print(tail)
print('ERR', r.stderr[-800:])
