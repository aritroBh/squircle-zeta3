import os
p = '<project>/research/08_theory/M4d_w5_nice'
print(open(os.path.join(p,'result.md')).read())
print("="*70)
print(open(os.path.join(p,'progress.log')).read())
