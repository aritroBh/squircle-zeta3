import os
work = '<project>/research/08_theory/WT2r_review'
os.makedirs(work, exist_ok=True)
base = '<project>/research/08_theory'
with open(base + '/WT2b_stabilization/result.md') as f:
    content = f.read()
print(content[:8000])
