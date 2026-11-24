import os, glob
# Current working directory
print("CWD:", os.getcwd())
print("\nContents of cwd:")
for p in sorted(glob.glob('*')):
    print(" ", p)
print("\nresearch/08_theory:")
for p in sorted(glob.glob('research/08_theory/*'))[:50]:
    print(" ", p)
