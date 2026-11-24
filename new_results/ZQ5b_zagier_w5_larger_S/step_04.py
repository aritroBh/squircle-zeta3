import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
# Read key files
with open('../ZQ5_zagier_weight5/w5lib.py') as f:
    print(f.read()[:6000])
