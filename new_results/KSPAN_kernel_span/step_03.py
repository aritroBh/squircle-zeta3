import os, json, sys
from pathlib import Path

cwd = Path.cwd()
print(f"Working in: {cwd}")
print("\nContents of current directory:")
for f in cwd.iterdir():
    print(f"  {f.name}")