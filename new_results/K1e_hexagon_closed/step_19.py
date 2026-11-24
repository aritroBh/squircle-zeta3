
import os, shutil, json, time
base='research/08_theory'
k1e=os.path.join(base,'K1e_hexagon_closed')
k1d=os.path.join(base,'K1d_engine_fixed')
k1c=os.path.join(base,'K1c_engine_predict')
for f in ['E3_fixed.pkl','Jtab_t.pkl','geometry.json','esum_engine.py']:
    shutil.copy(os.path.join(k1d,f), os.path.join(k1e,f))
for f in ['jloc.json','t3_b0_150.json','t3_b1_110.json']:
    shutil.copy(os.path.join(k1c,f), os.path.join(k1e,f))
print(sorted(os.listdir(k1e)))
with open(os.path.join(k1e,'progress.log'),'a') as fh:
    fh.write(f"# K1e progress log\n[{time.strftime('%H:%M:%S')}] step_01-05: exploration; read K1d result.md, K1c result.md, engine files; copied E3_fixed.pkl, Jtab_t.pkl, jloc.json, geometry.json, t3 files into work dir.\n")
