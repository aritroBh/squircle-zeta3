
import sympy as sp
import sym_G
G21, G12, syms = sym_G.build()
names = ['s0','s1','c0','c1','tau1','tau2','tau3','gam1','gam2','gam3','S0','Sc','lS0']
v = {k:syms[k] for k in names}

def monomials_in(expr, vars):
    # return set of exponent-tuples over vars appearing
    P = sp.Poly(sp.expand(expr), *vars)
    return set(P.monoms()), P

logvars = [v['tau1'],v['tau2'],v['tau3'],v['gam1'],v['gam2'],v['gam3'],v['lS0']]
for label,G in [('G21',G21),('G12',G12)]:
    mons,_ = monomials_in(G, logvars)
    nonzero = [mm for mm in mons if any(e!=0 for e in mm)]
    print(label, ": monomials in (tau1..3,gam1..3,lS0) with nonzero power:", nonzero)
    # So G depends on NONE of them?
    print("   depends only on s0,s1,c0,c1,S0,Sc :", len(nonzero)==0)
