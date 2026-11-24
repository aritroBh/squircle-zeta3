import os, pickle
import sympy as sp, mpmath as mp
WD = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
with open(f'{WD}/coef_acc_raw.pkl','rb') as f: coef_acc = pickle.load(f)
A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')

# Prop-8 local laws in terms of edge data (alpha=A0, beta=A1)
def Ze_laws(A0v,A1v,L2v,q0v,V02v):
    al,be = A0v,A1v
    u = -q0v/L2v; v = 1+q0v/L2v
    w = al*u+be*v; e = al*be*u*v; q = al**2*u**2+be**2*v**2
    Q = al**2+be**2; p = al*be; s = al+be
    h2 = V02v - q0v**2/L2v; rr = h2/L2v  # r = h^2/L^2
    h2overL = h2/sp.sqrt(L2v)
    Z5 = 360*al*be*h2overL*(20*w**3-20*e*w-32*q-48*e+15*w-2 + rr*(-15*Q*w+15*p*w+5*p*s+8*Q-12*p))
    Z23 = 120*al*be*h2overL*(-60*e*w+12*q+72*e-12*w+2 + rr*(-45*p*w+15*p*s-3*Q+18*p))
    return sp.simplify(Z5), sp.simplify(Z23)

# numeric general-r check at several random valid edge-data points
mp.mp.dps=50
import random
random.seed(1)
mkz5=((),0,(0,0,1)); mkz23=((),0,(1,1,0))
print("General-r check: 120*C_z5 == Z_e^(5), 120*C_z2z3 == Z_e^(2,3) ?")
worst5=mp.mpf(0); worst23=mp.mpf(0)
for trial in range(6):
    # random valid: pick B, kap, L2, q0, V02 with V02*L2>q0^2; A1=1-B, A0=1-B-B*kap
    Bv=sp.Rational(random.randint(1,4),random.randint(2,7))
    kv=sp.Rational(random.randint(-2,3),random.randint(2,5))
    L2v=sp.Rational(random.randint(2,9),1)
    q0v=sp.Rational(random.randint(-4,4),random.randint(1,3))
    V02v=sp.Rational(random.randint(2,9),1)
    A1v=1-Bv; A0v=1-Bv-Bv*kv
    if V02v*L2v<=q0v**2: continue
    rd={A0:A0v,A1:A1v,B:Bv,kap:kv,L2:L2v,q0:q0v,V02:V02v,lA0:sp.log(A0v),lA1:sp.log(A1v),lB:sp.log(Bv)}
    c5=120*coef_acc[mkz5].subs(rd); c23=120*coef_acc[mkz23].subs(rd)
    Z5,Z23=Ze_laws(A0v,A1v,L2v,q0v,V02v)
    d5=abs(mp.mpf(sp.N(c5-Z5,40))); d23=abs(mp.mpf(sp.N(c23-Z23,40)))
    worst5=max(worst5,d5); worst23=max(worst23,d23)
    print(f"  trial {trial}: |c5-Z5|={mp.nstr(d5,3)}, |c23-Z23|={mp.nstr(d23,3)}")
print(f"worst |c5-Z5|={mp.nstr(worst5,3)}, worst |c23-Z23|={mp.nstr(worst23,3)}")
