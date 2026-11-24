import mpmath as mp, pickle, sympy as sp
mp.mp.dps = 745
vals = pickle.load(open('basis745.pkl','rb'))
RED = pickle.load(open('reductions.pkl','rb'))
atoms = {k: mp.mpf(v) for k,v in pickle.load(open('oct5/atoms600.pkl','rb')).items()}
l2,lu,z2,z3 = sp.symbols('l2 lu z2 z3')
L2x1,L2x2,L2k,L2d,L2r2h,L2mr2h = sp.symbols('L2x1 L2x2 L2k L2d L2r2h L2mr2h')
L3x1,L3x2,L3k,L3d,L3r2h,L3mr2h = sp.symbols('L3x1 L3x2 L3k L3d L3r2h L3mr2h')
vmap = {l2:vals['l2'],lu:vals['lu'],z2:vals['z2'],z3:vals['z3'],
 L2x1:vals['Li2_x1'],L2x2:vals['Li2_x2'],L2k:vals['Li2_k1'],L2d:vals['Li2_d'],L2r2h:vals['Li2_r2h'],L2mr2h:vals['Li2_mr2h'],
 L3x1:vals['Li3_x1'],L3x2:vals['Li3_x2'],L3k:vals['Li3_k1'],L3d:vals['Li3_d'],L3r2h:vals['Li3_r2h'],L3mr2h:vals['Li3_mr2h']}
ev = lambda e: mp.mpf(sp.N(e.subs(vmap, simultaneous=True),740))
# stuffle: Li2(a)Li1(b) = Li21(a,b)+Li12(a,b)+Li3(ab)
# Li12_a = Li12(x1,-r2/2) = Li2(x1)*Li1(-r2/2) - Li3(x2) - Lam
Li1_mr2h = -(vals['lu']-vals['l2']/2)   # log(1+r2/2)?? careful: Li1(b) = -log(1-b); b=-r2/2 -> -log(1+r2/2)
Li1_mr2h = -mp.log(1+mp.sqrt(2)/2)
Li12_a_pred = vals['Li2_x1']*Li1_mr2h - vals['Li3_x2'] - ev(RED['Lam'])
print("Li12_a (stuffle+R3) =", mp.nstr(Li12_a_pred,25))
print("Li12_a (atoms600)   =", mp.nstr(atoms['Li12_a'],25))
print("|diff| =", mp.nstr(abs(Li12_a_pred-atoms['Li12_a']),3))
Li1_msqrt2 = -mp.log(1+mp.sqrt(2))
Li12_b_pred = vals['Li2_x2']*Li1_msqrt2 - vals['Li3_x1'] - ev(RED['swap'])
print("Li12_b (stuffle+R4) =", mp.nstr(Li12_b_pred,25))
print("Li12_b (atoms600)   =", mp.nstr(atoms['Li12_b'],25))
print("|diff| =", mp.nstr(abs(Li12_b_pred-atoms['Li12_b']),3))
# Li11 cross-check too: Li11(a,b) = Li1(a)Li1(b) - Li2(ab)... stuffle: Li1(a)Li1(b) = Li11(a,b)+Li11(b,a)+Li2(ab)
# with (a,b)=(x1,-r2/2): ab = x2; Li11(x1,-r2/2) = Li11_a. Known identity Li11(x,y)=Li2(x/(x-1)*...) skip.
# OCT4's Li3a cross-check: Li3(2*sqrt2-2) in my span?
from flint import fmpz_mat
vals2 = dict(vals); vals2['Li3a_oct4'] = mp.polylog(3, 2*(mp.sqrt(2)-1))
w3 = ['z3','Li3_mr2h','Li3_x1','Li3_x2','Li3_k1','Li3_d','Li3_r2h','z2l2','z2lu','l2^3','l2^2lu','l2lu^2','lu^3']
nb=len(w3); n=2*nb+1; C=mp.mpf(10)**600
vec=[vals2['Li3a_oct4']]
for b in w3: vec += [vals2[b], mp.sqrt(2)*vals2[b]]
rows=[]
for i in range(n):
    row=[0]*n+[int(mp.nint(C*vec[i]))]; row[i]=1; rows.append(row)
M=fmpz_mat(rows).lll()
v=[int(M[0,j]) for j in range(n)]
res=abs(sum(mp.mpf(v[j])*vec[j] for j in range(n)))
print("\nLi3(2k1) [OCT4's Li3a] in my span? maxc=%d resid=%s"%(max(abs(x) for x in v), mp.nstr(res,3)))
if res < mp.mpf('1e-300'):
    from fractions import Fraction as Fr
    for i,b in enumerate(w3):
        a,bb=v[1+2*i],v[2+2*i]
        if a or bb: print("   ", Fr(a,-v[0]), "+", Fr(bb,-v[0]), "r2  *", b)
