
import os
ROOT = 'research/08_theory/D13s_fix_recertify'
for dirpath, dirnames, filenames in os.walk(ROOT):
    depth = dirpath.replace(ROOT,'').count(os.sep)
    if depth > 1: continue
    print(dirpath)
    for f in sorted(filenames):
        print("   ", f, os.path.getsize(os.path.join(dirpath,f)))
