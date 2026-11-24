import glob, os
print("HEXC2 files:")
for p in sorted(glob.glob('../HEXC2_level6_certificate/*')):
    print("  ", os.path.basename(p), os.path.getsize(p))
print("\nCOR7 files:")
for p in sorted(glob.glob('../COR7_classical_corners/*')):
    print("  ", os.path.basename(p), os.path.getsize(p))
