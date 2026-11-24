
import shutil, os, time
chn1 = '<project>/research/08_theory/CHN1_chain_bodies'
for f in ['mzv_exact.py','engine.py','symengine.py','mzv_num.py','mzv_vals.py']:
    shutil.copy(f'{chn1}/{f}', f)
print("copied:", [f for f in os.listdir('.') if f.endswith('.py')])

# quick timing of DS tables at weights 2..10 (pure) and check W=11,12 flint timing
import mzv_exact as M
for W in range(2,11):
    t=time.time(); T=M.DSTable(W); print(f"W={W} words={T.n} nullity={T.nullity} {time.time()-t:.2f}s")
