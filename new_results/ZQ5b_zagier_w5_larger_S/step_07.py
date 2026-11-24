import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
with open('../ZQ5_zagier_weight5/pipe.py') as f:
    content = f.read()
print(f"pipe.py length: {len(content)}")
print(content[:8000])
