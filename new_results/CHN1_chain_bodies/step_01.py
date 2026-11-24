import os, shutil, subprocess, time
t0=time.time()
WORK='<project>/work_chn1'
os.makedirs(WORK, exist_ok=True)
os.chdir(WORK)
PILOT='<project>/research/02_ideas/lenses/experimental_pilots'
# list pilot dir
print(os.listdir(PILOT))
# copy the files we are told to reuse
for f in ['mzv_vals.py','chain4_generic_check.py','chain4_taylor.py','chain4_generic.py',
          'chain4_coeffs.txt']:
    src=os.path.join(PILOT,f)
    if os.path.exists(src):
        shutil.copy(src, os.path.join(WORK,f)); print('copied',f)
    else:
        print('MISSING',f)
print('PY ok:', os.path.exists('<project>/.venv/bin/python'))
print('t=%.1f'%(time.time()-t0))
