import os
# Find the work directory for this task
base = '<project>'
for root, dirs, files in os.walk(os.path.join(base, 'research')):
    if 'progress.log' in files:
        print(root)
        for f in sorted(files):
            print('   ', f)
