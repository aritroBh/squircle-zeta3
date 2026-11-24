import os
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
print(open('chain4_taylor.py').read())
print('#'*30)
print(open('chain4_generic.py').read()[:4000])
