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
import dec_engine as E
NAMED2 = [(x1, E.L2_1, E.L3_1), (x4, E.L2_4, E.L3_4), (1/phi, E.L2f, E.L3f)]
def canon_atoms(e):
    def replp(f):
        s_, arg = f.args[0], sp.simplify(f.args[1])
        for pt, s2, s3 in NAMED2:
            if sp.simplify(arg - pt) == 0:
                return s2 if s_ == 2 else s3
        raise ValueError(f"unknown polylog point: {arg}")
    e = e.replace(lambda f: isinstance(f, sp.Function) and f.func == sp.polylog and f.args[0] in (2,3) and not f.args[1].has(t), replp)
    # constant logs
    def repll(f):
        arg = sp.simplify(f.args[0])
        # candidates: phi powers, 2 powers
        for b in range(-6,7):
            if sp.simplify(arg - phi**b) == 0: return b*E.lf
        for a in range(-6,7):
            if sp.simplify(arg - 2**a) == 0: return a*E.l2
        raise ValueError(f"unknown log: {arg}")
    e = e.replace(lambda f: isinstance(f, sp.log) and not f.args[0].has(t), repll)
    return sp.expand(e)
res = {}
for name, F, Feng in [("D5",F5,E.D5), ("D6",F6,E.D6), ("D7",F7,E.D7), ("D8",F8,E.D8)]:
    val = sp.expand(F.subs(t,1) - F.subs(t,0))
    val = sp.expand(val.subs({sp.log(phi): E.lf, sp.log(2): E.l2}))
    vn = canon_atoms(val)
    dn = sp.expand(vn - Feng)
    print(name, "diff vs DEC1 engine:", dn)
    res[name] = vn
pickle.dump(res, open('D5678_proved.pkl','wb'))
