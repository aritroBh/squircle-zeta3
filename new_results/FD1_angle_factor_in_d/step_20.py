import pickle
fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
bs=pickle.load(open(fd1+'/beta_squircle.pkl','rb'))
ks=sorted(bs.keys())
print("squircle keys sample:", ks[:8], "...")
ds=sorted(set(k[0] for k in ks))
print("depths:", ds)
# print the d=5..8 top table entries
from fractions import Fraction as Fr
for d in [3,4,5,6,7,8]:
    row={k[1:]:v for k,v in bs.items() if k[0]==d}
    print(f"d={d}:", row)
