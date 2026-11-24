import os, subprocess
print("CWD:", os.getcwd())
print("--- cwd listing ---")
for f in sorted(os.listdir('.')):
    print(f)
print("--- parent listing ---")
for f in sorted(os.listdir('..')):
    print(f)
