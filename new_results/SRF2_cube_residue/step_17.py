import os
os.chdir("<project>/research/08_theory/SRF2_cube_residue")
for f in sorted(os.listdir(".")):
    if f.endswith(".py"):
        print(f, os.path.getsize(f))
