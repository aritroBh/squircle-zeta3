import mpmath as mp, pickle, sympy as sp
mp.mp.dps = 745
r2=mp.sqrt(2); k1=r2-1; x1=-k1
vals = pickle.load(open('basis745.pkl','rb'))
RED = pickle.load(open('reductions.pkl','rb'))
l2,lu,z2,z3 = sp.symbols('l2 lu z2 z3')
L2x1,L2x2,L2k,L2d,L2r2h,L2mr2h = sp.symbols('L2x1 L2x2 L2k L2d L2r2h L2mr2h')
L3x1,L3x2,L3k,L3d,L3r2h,L3mr2h = sp.symbols('L3x1 L3x2 L3k L3d L3r2h L3mr2h')
vmap = {l2:vals['l2'],lu:vals['lu'],z2:vals['z2'],z3:vals['z3'],
 L2x1:vals['Li2_x1'],L2x2:vals['Li2_x2'],L2k:vals['Li2_k1'],L2d:vals['Li2_d'],L2r2h:vals['Li2_r2h'],L2mr2h:vals['Li2_mr2h'],
 L3x1:vals['Li3_x1'],L3x2:vals['Li3_x2'],L3k:vals['Li3_k1'],L3d:vals['Li3_d'],L3r2h:vals['Li3_r2h'],L3mr2h:vals['Li3_mr2h']}
ev = lambda e: mp.mpf(sp.N(e.subs(vmap, simultaneous=True),740))
# Li12(a,b) = sum_{n>m>=1} a^n b^m/(n m^2), a=-r2/2, b=x1
def Li12ser(a,b):
    N = int(760/mp.log10(1/abs(a)))+60
    s=mp.mpf(0); ap=a; 
    for n in range(2,N+1):
        ap*=a
        inner=mp.mpf(0); bp=b
        for m in range(1,n):
            if m>1: bp*=b
            inner += bp/m**2
        s += ap*inner/n
    return s
import time; t0=time.time()
L12c = Li12ser(-r2/2, x1)
print("Li12(-r2/2,x1) series (%.1fs):"%(time.time()-t0), mp.nstr(L12c,25))
# stuffle: Li2(x1)*Li1(-r2/2) = Li21(x1,-r2/2) + Li12(-r2/2,x1) + Li3(x1*(-r2/2)=x2)
Li1_mr2h = -mp.log(1+r2/2)
lhs = vals['Li2_x1']*Li1_mr2h
rhs = ev(RED['Lam']) + L12c + vals['Li3_x2']
print("stuffle with R3: |Li2*Li1 - (Lam+Li12+Li3)| =", mp.nstr(abs(lhs-rhs),3))
# also with Lam direct series value:
rhs2 = vals['Lam'] + L12c + vals['Li3_x2']
print("stuffle with direct Lam:", mp.nstr(abs(lhs-rhs2),3))
