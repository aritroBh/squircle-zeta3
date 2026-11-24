import subprocess
# find the relevant scripts in W5C1
r = subprocess.run(['ls','../W5C1_certificate_scoping/*.py'], capture_output=True, text=True)
print(r.stdout[:200])
import glob
for f in sorted(glob.glob('../W5C1_certificate_scoping/*.py')):
    print(f, len(open(f).read()))