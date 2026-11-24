src = open('../HEXC2_level6_certificate/hexc2_common.py').read()
print(src[4300:6100])
print("=== HEXC1 files ===")
import glob, os
for p in sorted(glob.glob('../HEXC1*/*'))[:40]:
    print("  ", p)
