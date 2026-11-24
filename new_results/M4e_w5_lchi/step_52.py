import shutil, os
M4E = '<project>/research/08_theory/M4e_w5_lchi'
M4D = '<project>/research/08_theory/M4d_w5_nice'
for f in ['w5_verify.py','w5_cobracket.py']:
    shutil.copy(os.path.join(M4D,f), os.path.join(M4E,f))
print(open(M4E+'/w5_verify.py').read())
