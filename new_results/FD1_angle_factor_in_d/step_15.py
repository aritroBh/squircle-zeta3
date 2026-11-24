import shutil, os, pickle
a2='<project>/research/08_theory/ALPHA2_angle_leading_law'
fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
for f in ['corner_series.py','Fd_data.pkl','Fd_fits.pkl','beta_squircle.pkl','beta_tables.pkl','assemble.py',
          'cwfield.py','cwintegrand.py','cwred.py','cwbuild.py','lllkernel.py','hpl_tools.py']:
    shutil.copy(a2+'/'+f, fd1+'/'+f)
for f in os.listdir(a2):
    if f.startswith('kern18_') and f.endswith('.pkl'):
        shutil.copy(a2+'/'+f, fd1+'/'+f)
print(sorted(os.listdir(fd1)))
open(fd1+'/progress.log','a').write('step_10: copied corner_series, reducer machinery, kern18 pkls, Fd data into FD1\n')
