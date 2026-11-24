import os, subprocess
# First, figure out where we are and what's around
print("CWD:", os.getcwd())
print()
# List the parent research dir
for root in ['<project>']:
    print("=== ", root)
    try:
        for f in sorted(os.listdir(root)):
            print("  ", f)
    except Exception as e:
        print("ERR", e)
