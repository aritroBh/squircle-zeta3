import os, sys

# Check parent directory
parent = os.path.dirname(os.getcwd())
print(f"Parent directory: {parent}")
print("\nParent contents:")
for item in sorted(os.listdir(parent)):
    print(f"  {item}")