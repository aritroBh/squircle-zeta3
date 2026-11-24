import os, subprocess
base='<project>/research/08_theory'
print(os.listdir(base))
for d in ['M4c_w5_closed_form','M4e_w5_lchi','M4f_w5_extras','A3_W4_proof']:
    p=os.path.join(base,d)
    print('=====',d)
    print(sorted(os.listdir(p)))
