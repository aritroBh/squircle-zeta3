import os, glob
w=glob.glob('<project>/research/08_theory/WTS_weight_k_structure/**', recursive=True)
ks=sorted([f for f in w if 'kern' in f and f.endswith('.pkl')])
print(ks)
# also check LEAD1 folder for squircle kernels
l1='<project>/research/08_theory/LEAD1_leading_genfun'
print([f for f in os.listdir(l1) if f.endswith('.pkl') or 'kern' in f])
