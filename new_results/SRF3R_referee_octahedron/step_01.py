import os, subprocess
base = "research/08_theory/SRF3_octahedron_radius_proof"
cwd = os.getcwd()
print("CWD:", cwd)
for root, dirs, files in os.walk("research/08_theory"):
    level = root.count(os.sep)
    print("  "*level + os.path.basename(root) + "/")
    for f in files:
        p = os.path.join(root, f)
        print("  "*(level+1) + f, os.path.getsize(p))
