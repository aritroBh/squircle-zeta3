import os
base = "../SRF3_octahedron_radius_proof"
for f in sorted(os.listdir(base)):
    p = os.path.join(base, f)
    print(f, os.path.getsize(p))
print("---")
srf1 = "../SRF1_radii_3d"
print(os.listdir(srf1))
