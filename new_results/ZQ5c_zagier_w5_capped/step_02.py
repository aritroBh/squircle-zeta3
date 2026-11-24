import os
cwd = os.getcwd()
print("CWD:", cwd)
print("\nContents:")
for item in sorted(os.listdir(cwd)):
    p = os.path.join(cwd, item)
    if os.path.isdir(p):
        print(f"  {item}/")
    else:
        print(f"  {item}")
