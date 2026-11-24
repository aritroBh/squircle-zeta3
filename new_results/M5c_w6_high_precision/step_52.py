import os, glob
base='<project>/research/08_theory'
for d in ['M4e_w5_lchi','M4f_w5_close','M4i_w5_noL','M4h_doubled_basis_BRIEF.md']:
    p=os.path.join(base,d)
    print(d, os.path.exists(p))
print([d for d in os.listdir(base) if d.startswith('M4') and os.path.isdir(os.path.join(base,d))])
