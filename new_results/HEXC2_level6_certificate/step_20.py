import os
p='<project>/research/08_theory/COR4_M3alpha_classical'
print([f for f in os.listdir(p) if 'hpl' in f.lower()])
print(open(p+'/hpl_series.py').read())
