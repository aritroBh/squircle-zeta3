
import time
from fractions import Fraction as F
import mzv_exact as M
import engine5

t=time.time()
eng = engine5.build_engine(6)
print("engine built W<=6:", f"{time.time()-t:.1f}s")
for W in range(3,7):
    nw = len(eng['Areg'][W])
    print(f"W={W}: Areg words={nw}", "admissible:", all(w[0]>=2 and sum(w)==W for w in eng['Areg'][W]))
vw = engine5.volreg_words(eng, 6)
for W in range(3,7):
    print(f"W={W}: Vol_reg words={len(vw[W])}")
# reduced forms
t=time.time()
forms6, _ = engine5.reduced_forms(6)
print("reduced:", f"{time.time()-t:.1f}s")
import mzv_exact as MM
for W in range(3,7):
    T=engine5.table(W)
    print(f"W={W} basis={T.monomials}")
    for b,p in enumerate(forms6[W]):
        print(f"   coord{b}: {len(p)} terms", dict(list(p.items())[:4]))
