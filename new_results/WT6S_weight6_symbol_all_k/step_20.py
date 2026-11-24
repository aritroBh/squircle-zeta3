
import os
ROOT = '<project>'
T = os.path.join(ROOT, 'research/08_theory')
# find qsqrt2.py, w4_symbol.py, w4_right_coproduct.py
import glob
for name in ['qsqrt2.py', 'w4_symbol.py', 'w4_right_coproduct.py', 'w5_symbol.py']:
    hits = glob.glob(os.path.join(T, '*', name)) + glob.glob(os.path.join(ROOT, 'paper/squircle_note/scripts', name))
    print(name, "->", hits[:4])
print()
M5B = os.path.join(T, 'M5b_w6_closed_form')
print(sorted(os.listdir(M5B)))
