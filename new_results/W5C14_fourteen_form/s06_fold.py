"""Fold Delta (exact, word tuples) into the four Galois sectors with W5C2's exact fold code (read-only import)."""
import sys, pickle, time
sys.dont_write_bytecode=True
TH='<project>/research/08_theory'
sys.path.insert(0, TH+'/W5C2_exact_certificate')
from w5c2_fold import Ctx, fold_D, SECTORS
t0=time.monotonic()
ctx=Ctx(); assert ctx.check_vs_meta()
d=pickle.load(open('delta.pkl','rb'))
out={}
for comp in ('rat','sq'):
    dv={ctx.IDS[w]:a for w,a in d[comp].items()}
    out[comp]=dv
    for sec in SECTORS:
        out[(comp,sec)]=fold_D(dv,sec,ctx)
        print(comp,sec,'sector support',len(out[(comp,sec)]))
pickle.dump(out,open('delta_folds.pkl','wb'))
print(f'seconds={time.monotonic()-t0:.0f}')
