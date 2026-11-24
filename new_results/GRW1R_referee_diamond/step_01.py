import os
os.chdir('research/08_theory/GRW1R_referee_diamond') if os.path.exists('research/08_theory/GRW1R_referee_diamond') else None
for root, dirs, files in os.walk('.'):
    level = root.replace('.', '').count(os.sep)
    if level < 3:
        print('  '*level + os.path.basename(root) + '/')
        for f in files[:20]:
            print('  '*(level+1) + f)
