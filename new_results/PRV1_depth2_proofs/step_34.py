import os, glob

# Find the work directory. The task says "Work only in your work directory."
# Let's look around the research directory structure for my previous outputs.
base = '<project>'
for root, dirs, files in os.walk(base):
    # limit depth
    depth = root[len(base):].count(os.sep)
    if depth <= 3:
        if any(f.startswith('step_') or f in ('progress.log','result.md') for f in files):
            print(root)
            for f in sorted(files):
                print('   ', f)
