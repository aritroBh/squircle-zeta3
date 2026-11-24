import os
base = '<project>/research/08_theory'
print("Base:", base)
for item in sorted(os.listdir(base)):
    p = os.path.join(base, item)
    if os.path.isdir(p):
        print(f"  {item}/")
    else:
        print(f"  {item}")

# Check the current dir files
cwd = os.getcwd()
print("\n--- Current dir step files ---")
for f in ['step_01.py','step_02.py']:
    p = os.path.join(cwd, f)
    if os.path.exists(p):
        print(f"\n===== {f} =====")
        with open(p) as fh:
            print(fh.read())
