import sys, os, shutil, time
base = "<project>/research/08_theory/"
cwd = base+"DOD5_dodecagon_termwise/"
# copy libs for self-containment
for f in ["dod3_symbol.py","dod3_R.py","dod4_lib.py"]:
    shutil.copy(base+"DOD4_dodecagon_closed/"+f, cwd+f)
sys.path.insert(0, cwd)
t0=time.time()
from dod4_lib import *
print("import+init time: %.1fs"%(time.time()-t0))
print("BASIS3_ZS (8 Li3 args):")
for z in BASIS3_ZS: print("   ", z, "~", pfloat(z))
print("BASIS2_ZS:", BASIS2_ZS, [pfloat(z) for z in BASIS2_ZS])
t0=time.time()
C = consts(620)
print("consts time %.1fs"%(time.time()-t0))
t0=time.time()
vals = [Li3v(z) for z in BASIS3_ZS]
print("8 Li3v time %.1fs"%(time.time()-t0))
t0=time.time()
rem = mp.mpf('0.01382611404688009602330383739453034316291')
names, ivals = invisible10(C)
rel = mp.pslq([rem]+ivals, tol=mp.mpf('1e-565'), maxcoeff=10**12, maxsteps=200000)
print("pslq 620 time %.1fs, rel:"%(time.time()-t0), rel)
