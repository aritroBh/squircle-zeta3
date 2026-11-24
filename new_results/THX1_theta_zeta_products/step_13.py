
import os
for f in sorted(os.listdir('.')):
    if f.startswith('step_') or f == 'transcript.json':
        print('='*20, f, os.path.getsize(f))
        print(open(f).read()[:600])
