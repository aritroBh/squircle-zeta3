import os, shutil, json

SRC = 'research/08_theory/K4c_3d_edge_law/'
print("Files in K4c dir:")
for f in sorted(os.listdir(SRC)):
    print("  ", f, os.path.getsize(os.path.join(SRC,f)))

# Work directory
WORK = 'research/08_theory/K4d_zeta3_coeffs/'
os.makedirs(WORK, exist_ok=True)
print("\nWork dir:", WORK)
