import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
with open('../ZQ5_zagier_weight5/result.md') as f:
    content = f.read()
print(f"len {len(content)}")
print(content[:9000])
