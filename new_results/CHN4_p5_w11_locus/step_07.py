import sys
sys.path.insert(0, '../CHN3_p5_locus')
import mzv_exact as M
for W in range(2,13):
    print(f"W={W}: nullity(basis)={len(M.BASIS[W])}")
    for m in M.BASIS[W]:
        print("    ", m)
