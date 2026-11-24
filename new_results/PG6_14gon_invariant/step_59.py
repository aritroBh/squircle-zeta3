import os
print("CWD:", os.getcwd())
print("=== progress.log ===")
print(open('progress.log').read())
print("=== files ===")
for f in sorted(os.listdir('.')):
    if not f.startswith('step_'):
        print(" ", f, os.path.getsize(f))
