
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
g35 = make_f_eval(11,7)   # zeta3*zeta(3,5) coordinate
_ = g35((1,2,3,4,5))
def Lval(a):
    a1,a2,a3,a4,a5=a; return a2*a3*a4*(a1+a2+a3+a4+a5)
log=open('ctrl_z3z35.log','w')
def W(s): log.write(s+'\n'); log.flush()
# 1) divisibility by a2a3a4(Sum a)?
tests={'a2=0':(1,0,2,3,1),'a3=0':(1,2,0,3,1),'a4=0':(1,2,3,0,1),'Sum=0':(1,2,3,4,-10)}
for nm,p in tests.items():
    W(f"{nm}: g35={g35(p)}")
# 2) interpolate deg-7 quotient h = g35/(a2a3a4 Sum a) using the SAME 170 det points as Q11
d=pickle.load(open('interp_points.pkl','rb')); sel=[tuple(x) for x in d['sel']]; orbits=d['orbits']; mons7=d['mons7']
orbit_of={}
for i,o in enumerate(orbits):
    for e in o: orbit_of[e]=i
NORB=len(orbits)
def orbit_vals_int(a):
    vals=[0]*NORB
    for e in mons7:
        pw=1
        for x,ei in zip(a,e):
            if ei: pw*=x**ei
        vals[orbit_of[e]]+=pw
    return vals
from flint import fmpq_mat, fmpq
Mdata=[]; rhs=[]
for a in sel:
    v=g35(a)
    Mdata += [fmpq(x) for x in orbit_vals_int(a)]
    rhs.append(fmpq(v.numerator,v.denominator)/Lval(a))
Mm=fmpq_mat(NORB,NORB,Mdata); bm=fmpq_mat(NORB,1,rhs)
sol=Mm.solve(bm)
coefs=[F(int(sol[i,0].p),int(sol[i,0].q)) for i in range(NORB)]
W(f"interpolated h: nonzero orbit coefs={sum(1 for c in coefs if c!=0)}, denom lcm check")
# build explicit quotient poly H7
H7={}
for i,o in enumerate(orbits):
    c=coefs[i]
    if c==0: continue
    for e in o: H7[e]=H7.get(e,F(0))+c
def peval(poly,a):
    a=[F(x) for x in a]; tot=F(0)
    for e,c in poly.items():
        pw=F(1)
        for x,ei in zip(a,e):
            if ei: pw*=x**ei
        tot+=c*pw
    return tot
# verify g35 = L*H7 on 20 fresh points
import random
random.seed(77)
bad=0
for _ in range(20):
    a=tuple(random.randint(1,9) for _ in range(5))
    if g35(a)!=Lval(a)*peval(H7,a): bad+=1
W(f"verify g35 = a2a3a4(Sum a)*H7 on 20 fresh pts: bad={bad}")
pickle.dump(H7, open('H7_z3z35.pkl','wb'))
# 3) is H7 divisible by C4?  (H7 deg7, C4 deg4 -> quotient deg3) exact division
import sympy as sp
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5'); A=[a1,a2,a3,a4,a5]
H7sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod(A[i]**e[i] for i in range(5)) for e,c in H7.items())
C4=sp.sympify(open('C4.txt').read())
q,r=sp.div(H7sp, C4, a1,a2,a3,a4,a5)
W(f"H7 / C4: remainder == 0 ? {sp.expand(r)==0}")
W(f"H7 degree={sp.Poly(H7sp,*A).total_degree()}")
if sp.expand(r)==0:
    W(f"quotient (deg3) terms={len(sp.Add.make_args(sp.expand(q)))}")
    pickle.dump(str(sp.expand(q)), open('H7_div_C4_quotient.pkl','wb'))
W("DONE")
log.close()
print("done")
