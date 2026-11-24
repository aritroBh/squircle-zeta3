import sympy as sp, pickle, os, sys
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd); sys.path.insert(0, '<project>/research/08_theory/DEC1_decagon')
t = sp.Symbol('t'); r5 = sp.sqrt(5); phi = (1+r5)/2
x1 = -1/phi; x4 = 1/phi**2; kap1 = 1/phi; kap4 = -1/phi**2
u1 = 1+kap1*t; u4 = 1+kap4*t
w1 = sp.log(u1); w4 = sp.log(u4)
lf = sp.log(phi); l2 = sp.log(2)
Li2 = lambda z: sp.polylog(2,z); Li3 = lambda z: sp.polylog(3,z)
F5 = phi*(lf*w1 - Li2(x4*(1+t/phi)))
F6 = -phi**2*(2*lf*w4 - Li2((1/phi)*(1+kap4*t)))
s7 = 1+kap4*t
G7 = 4*lf**2*sp.log(s7) - 4*lf*Li2(s7/phi) + sp.log(s7/phi)*sp.log(1-s7/phi)**2 \
     + 2*sp.log(1-s7/phi)*Li2(1-s7/phi) - 2*Li3(1-s7/phi)
F7 = -phi**2*G7
s8 = 1+t/phi
H8 = lf**2*sp.log(s8) - 2*lf*Li2(x4*s8) + sp.log(x4*s8)*sp.log(1-x4*s8)**2 \
     + 2*sp.log(1-x4*s8)*Li2(1-x4*s8) - 2*Li3(1-x4*s8)
F8 = phi*H8
# exact endpoints (all arguments inside (0,1) -> direct substitution is exact)
def endpoints(F):
    return sp.expand(F.subs(t,1)), sp.expand(F.subs(t,0))
import dec_engine as E
# named symbols in engine
def to_named(e):
    return e.subs({sp.log(phi): E.lf,
                   sp.polylog(2,x4): E.L2_4, sp.polylog(3,x4): E.L3_4,
                   sp.polylog(2,1/phi): E.L2f, sp.polylog(3,1/phi): E.L3f}, simultaneous=True)
res = {}
for name, F, Feng in [("D5",F5,E.D5), ("D6",F6,E.D6), ("D7",F7,E.D7), ("D8",F8,E.D8)]:
    F1, F0 = endpoints(F)
    val = sp.expand(F1 - F0)
    val = sp.expand(val.subs({sp.log(phi): E.lf}))
    # engine D in named atoms
    d = sp.expand(val - sp.expand(Feng.subs({E.lf: E.lf})))
    # replace remaining polylogs by engine symbols both sides: do full substitution
    vn = to_named(val)
    dn = sp.expand(vn - Feng)
    print(name, "endpoint value computed; symbolic diff vs DEC1 engine:", dn)
    res[name] = vn
pickle.dump(res, open('D5678_proved.pkl','wb'))
