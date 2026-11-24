import os
print("CWD:", os.getcwd())
print("Contents:")
for item in sorted(os.listdir(".")):
    p = os.path.join(".", item)
    print(f"  {'DIR ' if os.path.isdir(p) else 'FILE'} {item}")
