import glob, os
for f in sorted(glob.glob('step_*.py')):
    head = open(f).read(200).split('\n')[0]
    print(f"{f:16s} {os.path.getsize(f):6d}  {head[:100]}")
