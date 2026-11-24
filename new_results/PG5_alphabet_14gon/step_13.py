import os, shutil
# set up workspace
os.makedirs('.', exist_ok=True)
shutil.copy('../PGEN1R_referee_alphabet/symbol_engine_r.py', 'symbol_engine_r.py')
open('progress.log','a').write('PG5 start: copied symbol_engine_r.py from PGEN1R\n')

# find the paper tex
import glob
for pat in ['../../09_paper_v2/**/*.tex', '../../../09_paper_v2/**/*.tex', '../../09_paper*', '../../**/main.tex']:
    hits = glob.glob(pat, recursive=True)
    if hits:
        print(pat, '->', hits[:20])
