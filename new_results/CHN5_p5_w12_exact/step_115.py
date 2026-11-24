import os, sys, time, pickle
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
import point_engine as PE, mzv_exact as M
from fractions import Fraction as F
import random
Q=pickle.load(open('Q8_quotients.pkl','rb'))
def peval(Qd,a):
    a=[F(x) for x in a]; tot=F(0)
    for e,c in Qd.items():
        pw=F(1)
        for x,ei in zip(a,e):
            if ei: pw*=x**ei
        tot+=c*pw
    return tot
def Lval(a):
    a1,a2,a3,a4,a5=a; return F(a2)*a3*a4*(a1+a2+a3+a4+a5)
T=PE.table(12)
def engine_coords(pt):
    a=[F(x) for x in pt]; g=[F(0)]*12
    for p in (a,a[::-1]):
        vw=PE.volreg_words_point(p,12)
        for w,pp in vw[12].items():
            rw=T.reduce_word(w)
            for b in range(12):
                if rw[b]: g[b]+=pp*rw[b]
    return g
random.seed(777)
t0=time.time()
for trial in range(3):
    pt=tuple(random.randint(1,9) for _ in range(5))
    g=engine_coords(pt)
    L=Lval(pt)
    predA=L*peval(Q['A'],pt); predB=L*peval(Q['B'],pt)
    predC=L*peval(Q['CLASS'],pt); predZ=L*peval(Q['Z'],pt)
    okA=(g[7]==predA); okB=(g[8]==predB)
    okC=(g[9]+F(14,9)*g[10]==predC); okZ=(g[11]==predZ)
    print(f"pt {pt}: A match={okA}, B match={okB}, CLASS match={okC}, Z match={okZ}  [{time.time()-t0:.0f}s]")
