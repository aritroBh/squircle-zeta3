print(open('progress.log').read())
import os
files = sorted(os.listdir('.'))
print("deliverables:", [f for f in files if f.endswith('.md') or f.endswith('.py') and 'verify' in f])
print("verify_grw1r.py size:", os.path.getsize('verify_grw1r.py'))
