
import sys, time, pickle
sys.path.insert(0,'.')
import point_engine as PE, mzv_exact as M
from fractions import Fraction as F
def make_f_eval(W, idx):
    T = PE.table(W)
    def f(a):
        a=[F(x) for x in a]; acc=F(0)
        for pt in (a, a[::-1]):
            vw = PE.volreg_words_point(pt, W)
            for w,pp in vw[W].items():
                c=T.reduce_word(w)[idx]
                if c: acc+=pp*c
        return acc
    return f
f11 = make_f_eval(11,8)
_ = f11((1,2,3,4,5))
dd=pickle.load(open('det11_points.pkl','rb'))
sel=[tuple(x) for x in dd['sel']]
orbits11=dd['orbits11']; mons11=dd['mons11']
Q11=pickle.load(open('Q11_poly.pkl','rb'))
def Lval(a):
    a1,a2,a3,a4,a5=a; return a2*a3*a4*(a1+a2+a3+a4+a5)
def peval(poly,a):
    a=[F(x) for x in a]; tot=F(0)
    for e,c in poly.items():
        pw=F(1)
        for x,ei in zip(a,e):
            if ei: pw*=x**ei
        tot+=c*pw
    return tot
with open('proof_div.log','w') as lg:
    bad=0; t0=time.time()
    for i,a in enumerate(sel):
        D=f11(a)-Lval(a)*peval(Q11,a)
        if D!=0: bad+=1; lg.write(f"NONZERO at {a}: {D}\n")
        if (i+1)%50==0: lg.write(f"{i+1}/{len(sel)} checked bad={bad} {time.time()-t0:.0f}s\n"); lg.flush()
    lg.write(f"DONE: {len(sel)} points, {bad} nonzeros, {time.time()-t0:.0f}s\n")
    lg.write("VERDICT: "+("f11 = a2a3a4(Sum a)Q11 PROVED (D==0 on rank-693 det set)" if bad==0 else "FAIL")+"\n")
print("proof done")
