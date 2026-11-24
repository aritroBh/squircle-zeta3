import os, shutil, sys
from pathlib import Path

cwd = Path.cwd()
wt2_path = Path("../WT2_weight2_all_k")

# Copy binteg pickles for k=2..14 to current directory
for k in range(2, 15):
    src = wt2_path / f"binteg_{k}.pkl"
    dst = cwd / f"binteg_{k}.pkl"
    if src.exists():
        shutil.copy2(src, dst)
        print(f"Copied {src.name} to {dst}")
    else:
        print(f"Warning: {src} doesn't exist")

print(f"\nNow in {cwd}:")
for f in sorted(cwd.iterdir()):
    if f.name.startswith('binteg_'):
        print(f"  {f.name}")