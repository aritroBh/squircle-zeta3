
import os, sys, time, pickle
os.chdir('<project>/research/08_theory/CHN3_p5_locus')
sys.path.insert(0, '.')
import engine5, mzv_exact as M
t0=time.time()
eng = engine5.build_engine(12)
print(f"build_engine(12): {time.time()-t0:.1f}s", flush=True)
pickle.dump(eng, open('eng12.pkl','wb'))
t0=time.time()
vw = engine5.volreg_words(eng, 12)
print(f"volreg_words: {time.time()-t0:.1f}s; word counts W11={len(vw[11])} W12={len(vw[12])}", flush=True)
pickle.dump(vw, open('vw12.pkl','wb'))
for W in (11,12):
    t0=time.time()
    T = engine5.table(W)
    print(f"table({W}): {time.time()-t0:.1f}s nullity={T.nullity}", flush=True)
print("DONE", flush=True)
