
import os, sys, time, pickle
os.chdir('<project>/research/08_theory/CHN3_p5_locus')
sys.path.insert(0, '.')
import engine5
t0=time.time()
eng = engine5.build_engine(12)
print(f"build_engine(12): {time.time()-t0:.1f}s", flush=True)
pickle.dump(eng, open('eng12.pkl','wb'))
t0=time.time()
vw = engine5.volreg_words(eng, 12)
print(f"volreg_words: {time.time()-t0:.1f}s", flush=True)
pickle.dump(vw, open('vw12.pkl','wb'))
forms={}
for W in (11,12):
    t0=time.time()
    T=engine5.table(W); d=T.nullity
    acc=[{} for _ in range(d)]
    for w,pp in vw[W].items():
        rw=T.reduce_word(w)
        for b in range(d):
            if rw[b]: acc[b]=engine5.padd(acc[b], engine5.pscale(pp, rw[b]))
    forms[W]=acc
    print(f"reduce W={W}: {time.time()-t0:.1f}s", flush=True)
pickle.dump(forms, open('forms1112.pkl','wb'))
print("DONE", flush=True)
