
import os, sys, time, pickle
os.chdir('<project>/research/08_theory/CHN3_p5_locus')
sys.path.insert(0, '.')
import engine5
t0=time.time()
forms, eng = engine5.reduced_forms(12)
print(f"built forms W<=12 in {time.time()-t0:.1f}s", flush=True)
pickle.dump(forms, open('forms12.pkl','wb'))
for W in sorted(forms):
    print(W, [ (b, len(forms[W][b])) for b in range(len(forms[W])) if forms[W][b] ], flush=True)
print("DONE", flush=True)
