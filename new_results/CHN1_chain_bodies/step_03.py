import os
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
for f in ['mzv_vals.py','chain4_generic_check.py']:
    print('='*20, f, '='*20)
    print(open(f).read())
