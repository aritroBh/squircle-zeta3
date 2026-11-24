
import subprocess, os
print("CWD:", os.getcwd())
print("Files in CWD:")
for f in sorted(os.listdir('.')):
    print("  ", f)
print()
# Look for the work directory structure
base = '<project>'
if os.path.exists(base):
    print("Contents of", base)
    for f in sorted(os.listdir(base)):
        print("  ", f)
