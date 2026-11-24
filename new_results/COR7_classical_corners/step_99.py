import subprocess, sys
r = subprocess.run([sys.executable, 'verify_cor7.py'], capture_output=True, text=True, timeout=300)
print(r.stdout.splitlines()[-2])
print("returncode:", r.returncode)
import os
print(sorted(f for f in os.listdir('.') if not f.startswith('step_') and f != 'transcript.json'))
