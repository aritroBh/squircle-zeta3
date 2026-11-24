import shutil, os
base = '<project>/research/08_theory/'
pg5 = base+'PG5_alphabet_14gon/'
cwd = os.getcwd()
for fn in ['ftfield.py','geom14.pkl','E3n14.pkl','letters14.pkl','gold14.pkl']:
    shutil.copy(pg5+fn, os.path.join(cwd,fn))
print("copied:", os.listdir(cwd))
# Write initial progress.log and summary.md
open('progress.log','a').write('PG6 start: resolved to find P3 discrepancy; copied ftfield/geom14/E3n14/letters14/gold14 from PG5\n')
print(open('progress.log').read())
