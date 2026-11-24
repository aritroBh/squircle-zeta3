import os
# Post-interruption state check: confirm all deliverables are in place
files = sorted(os.listdir('.'))
need = ['result.md', 'verify_zlad1r.py', 'progress.log', 'muR_exact_N40_K8.json', 'muR_exact_N10_K9.json']
for f in need:
    print(f, os.path.exists(f), os.path.getsize(f) if os.path.exists(f) else '-')
print()
print('--- tail of progress.log ---')
print('\n'.join(open('progress.log').read().strip().split('\n')[-6:]))
