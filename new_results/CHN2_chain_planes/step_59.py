
import time, pickle, engine5, kimilib
from fractions import Fraction as F
import engine as E4

t=time.time(); eng10 = engine5.build_engine(10); print(f"build_engine(10): {time.time()-t:.1f}s", flush=True)
vw10 = engine5.volreg_words(eng10, 10)
print("words per W:", [len(vw10[W]) for W in range(3,11)])
forms10 = pickle.load(open('forms8_W8.pkl','rb'))
for W in range(9,11):
    t=time.time()
    T=engine5.table(W); d=T.nullity
    acc=[{} for _ in range(d)]
    for w,pp in vw10[W].items():
        rw=T.reduce_word(w)
        for b in range(d):
            if rw[b]: acc[b]=engine5.padd(acc[b], engine5.pscale(pp, rw[b]))
    forms10[W]=acc
    print(f"W={W} reduced in {time.time()-t:.1f}s", flush=True)
pickle.dump(forms10, open('forms10.pkl','wb'))
print("W=9 basis:", engine5.table(9).monomials)
print("W=10 basis:", engine5.table(10).monomials)
