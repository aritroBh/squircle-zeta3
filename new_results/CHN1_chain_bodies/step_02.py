import os, shutil, time
t0=time.time()
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
PILOT='<project>/research/02_ideas/lenses/experimental_pilots'
print(sorted(os.listdir(PILOT)))
for f in ['mzv_vals.py','chain4_generic_check.py','chain4_taylor.py','chain4_generic.py','chain4_coeffs.txt']:
    src=os.path.join(PILOT,f)
    if os.path.exists(src):
        shutil.copy(src, os.path.join(WORK,f)); print('copied',f)
    else:
        print('MISSING',f)
print('cwd writable test:', open('._w','w') and os.remove('._w') or 'ok')
print('t=%.1f'%(time.time()-t0))
