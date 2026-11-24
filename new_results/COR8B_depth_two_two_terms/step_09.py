
import subprocess
r = subprocess.run(['ps','aux'],capture_output=True,text=True)
lines = [l for l in r.stdout.splitlines() if 'python' in l.lower() or 'nice' in l.lower()]
print('\n'.join(lines))
