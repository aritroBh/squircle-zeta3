import os

# First, explore the directory structure
base = '<project>'
print("=== Current working directory ===")
print(os.getcwd())
print()
print("=== Contents of base ===")
for d in sorted(os.listdir(base)):
    print(d)
