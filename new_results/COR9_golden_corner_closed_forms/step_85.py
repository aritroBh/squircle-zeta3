WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
src = open(WD+'verify_cor9.py').read()
for line in src.splitlines():
    if 'B4 W3' in line or 'B5 W4' in line: print(repr(line))
