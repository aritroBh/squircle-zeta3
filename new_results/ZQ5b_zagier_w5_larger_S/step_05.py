import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
with open('../ZQ5_zagier_weight5/w5lib.py') as f:
    content = f.read()
print(f"total length: {len(content)}")
print(content[6000:14000])
