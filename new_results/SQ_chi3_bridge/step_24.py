import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C

# Check regularity of individual words in DIFF
for w in sorted(DIFF.keys()):
    reg_w = C.reg(w)
    if not (len(reg_w) == 1 and list(reg_w.keys())[0] == w):
        print(f"Word {w} regularizes to {reg_w}")

print("Done checking")