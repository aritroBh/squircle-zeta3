import os
print(sorted(os.listdir('.')))
for f in sorted(os.listdir('.')):
    if f.endswith('.py'):
        print('---', f, 'first line:', open(f).readline().strip()[:80])
