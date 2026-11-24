import subprocess, os
A1='<project>/research/08_theory/ALPHA1_angle_family_across_k'
A2='<project>/research/08_theory/ALPHA2_angle_leading_law'
for f in ['cwfield.py','cwintegrand.py','cwred.py','cwbuild.py','lllkernel.py','hpl_tools.py','corner_series.py','beta_tables.pkl','beta_squircle.pkl','kern_c12.pkl','kern_c35.pkl','kern_c13.pkl','kern_c45.pkl']:
    subprocess.run(['cp', f'{A1}/{f}', f'{A2}/{f}'])
print(subprocess.run(['ls','-la',A2], capture_output=True, text=True).stdout)
with open(f'{A2}/progress.log','a') as fh: fh.write('step_05: copied ALPHA1 reducer+data into ALPHA2\n')
