base = "<project>/research/08_theory/"
import glob, os
# find symeng.py
for f in glob.glob(base+"E5c_nonlocal_z5/*"):
    print(f, os.path.getsize(f) if os.path.isfile(f) else "<dir>")
