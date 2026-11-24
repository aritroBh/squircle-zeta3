import glob
for f in sorted(glob.glob('step_*.py')):
    c = open(f).read()
    if 'lA1' in c or 'unweighted' in c or 'logpq' in c:
        print("="*25, f)
        for i,line in enumerate(c.split('\n')):
            if 'lA1' in line or 'unweighted' in line or 'logpq' in line:
                print(f"  {i:3}: {line}")
