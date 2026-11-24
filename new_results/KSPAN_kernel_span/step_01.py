import os, sys, json, pickle, math, fractions
from pathlib import Path

# Check current directory and contents
print("Current directory:", os.getcwd())
print("\nDirectory contents:")
for item in os.listdir('.'):
    print(f"  {item}")