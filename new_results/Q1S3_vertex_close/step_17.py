import os
for f in sorted(os.listdir('.')):
    if f.endswith('.py'):
        print("="*70); print(f)
        print(open(f).read()[:600])
