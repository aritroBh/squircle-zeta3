import pickle, mpmath as mp, sympy as sp
from flint import fmpz_mat
from fractions import Fraction as Fr
mp.mp.dps = 745
vals = pickle.load(open('basis745.pkl','rb'))
r2m = mp.sqrt(2)
w3pure = ['z3','Li3_mr2h','Li3_x1','Li3_x2','Li3_k1','Li3_d','Li3_r2h','z2l2','z2lu','l2^3','l2^2lu','l2lu^2','lu^3',
      'Li2_mr2h*l2','Li2_mr2h*lu','Li2_x1*l2','Li2_x1*lu','Li2_x2*l2','Li2_x2*lu','Li2_k1*l2','Li2_k1*lu',
      'Li2_d*l2','Li2_d*lu','Li2_r2h*l2','Li2_r2h*lu']
def lll_first(target, basis_names, dps_scale=0.78):
    nb=len(basis_names); n=2*nb+1
    C=mp.mpf(10)**int(mp.mp.dps*dps_scale)
    vec=[vals[target]]
    for b in basis_names: vec += [vals[b], r2m*vals[b]]
    rows=[]
    for i in range(n):
        row=[0]*n+[int(mp.nint(C*vec[i]))]; row[i]=1; rows.append(row)
    M=fmpz_mat(rows).lll()
    cands=[]
    for i in range(n):
        v=[int(M[i,j]) for j in range(n)]
        if v[0]==0: continue
        res=abs(sum(mp.mpf(v[j])*vec[j] for j in range(n))); maxc=max(abs(x) for x in v)
        cands.append((maxc,res,v))
    cands.sort(key=lambda z:z[0]); return cands[0]
r2 = sp.sqrt(2)
sym = {'z3':sp.symbols('z3'),'Li3_mr2h':sp.symbols('L3mr2h'),'Li3_x1':sp.symbols('L3x1'),'Li3_x2':sp.symbols('L3x2'),
 'Li3_k1':sp.symbols('L3k'),'Li3_d':sp.symbols('L3d'),'Li3_r2h':sp.symbols('L3r2h'),
 'z2l2':sp.symbols('z2')*sp.symbols('l2'),'z2lu':sp.symbols('z2')*sp.symbols('lu'),
 'l2^3':sp.symbols('l2')**3,'l2^2lu':sp.symbols('l2')**2*sp.symbols('lu'),'l2lu^2':sp.symbols('l2')*sp.symbols('lu')**2,
 'lu^3':sp.symbols('lu')**3,
 'Li2_mr2h*l2':sp.symbols('L2mr2h')*sp.symbols('l2'),'Li2_mr2h*lu':sp.symbols('L2mr2h')*sp.symbols('lu'),
 'Li2_x1*l2':sp.symbols('L2x1')*sp.symbols('l2'),'Li2_x1*lu':sp.symbols('L2x1')*sp.symbols('lu'),
 'Li2_x2*l2':sp.symbols('L2x2')*sp.symbols('l2'),'Li2_x2*lu':sp.symbols('L2x2')*sp.symbols('lu'),
 'Li2_k1*l2':sp.symbols('L2k')*sp.symbols('l2'),'Li2_k1*lu':sp.symbols('L2k')*sp.symbols('lu'),
 'Li2_d*l2':sp.symbols('L2d')*sp.symbols('l2'),'Li2_d*lu':sp.symbols('L2d')*sp.symbols('lu'),
 'Li2_r2h*l2':sp.symbols('L2r2h')*sp.symbols('l2'),'Li2_r2h*lu':sp.symbols('L2r2h')*sp.symbols('lu')}
maxc,res,v = lll_first('swap', w3pure)
print("maxc",maxc,"resid",mp.nstr(res,3))
expr = sum((sp.Rational(v[1+2*i],-v[0]) + sp.Rational(v[2+2*i],-v[0])*r2)*sym[b] for i,b in enumerate(w3pure))
expr = sp.expand(expr)
print("rebuilt R4 terms:")
for term in sp.Add.make_args(expr): print("   ", term)
vmap = {sp.symbols('l2'):vals['l2'],sp.symbols('lu'):vals['lu'],sp.symbols('z2'):vals['z2'],sp.symbols('z3'):vals['z3'],
 sp.symbols('L2x1'):vals['Li2_x1'],sp.symbols('L2x2'):vals['Li2_x2'],sp.symbols('L2k'):vals['Li2_k1'],sp.symbols('L2d'):vals['Li2_d'],
 sp.symbols('L2r2h'):vals['Li2_r2h'],sp.symbols('L2mr2h'):vals['Li2_mr2h'],
 sp.symbols('L3x1'):vals['Li3_x1'],sp.symbols('L3x2'):vals['Li3_x2'],sp.symbols('L3k'):vals['Li3_k1'],sp.symbols('L3d'):vals['Li3_d'],
 sp.symbols('L3r2h'):vals['Li3_r2h'],sp.symbols('L3mr2h'):vals['Li3_mr2h']}
rhs = mp.mpf(sp.N(expr.subs(vmap, simultaneous=True),740))
print("|swap - R4_rebuilt| =", mp.nstr(abs(vals['swap']-rhs),3))
RED = pickle.load(open('reductions.pkl','rb'))
RED['swap'] = expr
pickle.dump(RED, open('reductions.pkl','wb'))
print("reductions.pkl updated")
