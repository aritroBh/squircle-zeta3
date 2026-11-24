
import time, engine5, kimilib
from fractions import Fraction as F
import engine as E4

t=time.time(); eng8 = engine5.build_engine(8); print(f"build_engine(8): {time.time()-t:.1f}s")
t=time.time()
vw8 = engine5.volreg_words(eng8, 8)
print(f"volreg: {time.time()-t:.1f}s; words per W: {[len(vw8[W]) for W in range(3,9)]}")
t=time.time()
forms8 = {}
for W in range(3,9):
    T=engine5.table(W); d=T.nullity
    acc=[{} for _ in range(d)]
    for w,pp in vw8[W].items():
        rw=T.reduce_word(w)
        for b in range(d):
            if rw[b]: acc[b]=engine5.padd(acc[b], engine5.pscale(pp, rw[b]))
    forms8[W]=acc
    print(f"W={W} reduced, {time.time()-t:.1f}s total")
import pickle
pickle.dump(forms8, open('forms8_W8.pkl','wb'))
# quick controls at W=7,8: u1=0 vs P4
def full_at(a, W):
    T=engine5.table(W); d=T.nullity
    reg=[kimilib.peval5(forms8[W][b], a) for b in range(d)]
    m0=kimilib.m0_coords(a, W)[W]
    return [reg[b]+m0[b] for b in range(d)]
ok=True
for pt in [(1,1,1,1,1),(2,1,1,3,1),(1,2,3,2,1)]:
    a1,a2,a3,a4,a5 = pt
    for W in (7,8):
        if full_at((0,a2,a3,a4,a5), W) != E4.reduced_coeffs((a2,a3,a4,a5), W)[W]:
            ok=False; print("u1=0 FAIL", pt, W)
        if full_at((a1,a2,a3,a4,0), W) != E4.reduced_coeffs((a1,a2,a3,a4), W)[W]:
            ok=False; print("u5=0 FAIL", pt, W)
print("W=7,8 degeneration controls:", "PASS" if ok else "FAIL")
print("W=8 basis:", engine5.table(8).monomials)
print("z(3,5) coord (b=3) terms:", len(forms8[8][3]))
