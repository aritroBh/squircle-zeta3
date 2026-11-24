import pickle, mpmath as mp, sympy as sp
mp.mp.dps = 780
r2=mp.sqrt(2); k1=r2-1; k2=r2/2-1; x1=-k1; x2=-k2
def Li21N(a,b,N):
    s=mp.mpf(0); ap=a*a; inner=b; bp=b
    for n in range(2,N+1):
        if n>2: ap*=a; bp*=b; inner+=bp/(n-1)
        s += ap*inner/n**2
    return s
# swap: a=x2, b=-r2 -> effective ratio |x2*r2|=k1: need N ~ 780/log10(1/k1) ~ 2037
import time; t0=time.time()
swap_hi = Li21N(x2,-r2,2200)
print("swap recomputed %.1fs"%(time.time()-t0), mp.nstr(swap_hi,30))
vals = pickle.load(open('basis745.pkl','rb'))
print("diff vs basis745 stored:", mp.nstr(abs(swap_hi-vals['swap']),3))
vals['swap'] = swap_hi
pickle.dump(vals, open('basis745.pkl','wb'))
# also recompute Lam at same rigor (effective ratio x1*r2/2 = (2-r2)/2=0.293; N~1500 enough) for safety
Lam_hi = Li21N(x1,-r2/2,1700)
print("Lam diff vs stored:", mp.nstr(abs(Lam_hi-vals['Lam']),3))
vals['Lam']=Lam_hi; pickle.dump(vals, open('basis745.pkl','wb'))
# save swap relation candidates at higher confidence: rerun swap LLL
from flint import fmpz_mat
from fractions import Fraction as Fr
def lll_rel(target, basis_names, dps_scale=0.78):
    nb=len(basis_names); n=2*nb+1
    C=mp.mpf(10)**int(mp.mp.dps*dps_scale)
    vec=[vals[target]]
    for b in basis_names: vec += [vals[b], r2*vals[b]]
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
    cands.sort(key=lambda z:z[0]); return cands
w3pure = ['z3','Li3_mr2h','Li3_x1','Li3_x2','Li3_k1','Li3_d','Li3_r2h','z2l2','z2lu','l2^3','l2^2lu','l2lu^2','lu^3',
      'Li2_mr2h*l2','Li2_mr2h*lu','Li2_x1*l2','Li2_x1*lu','Li2_x2*l2','Li2_x2*lu','Li2_k1*l2','Li2_k1*lu',
      'Li2_d*l2','Li2_d*lu','Li2_r2h*l2','Li2_r2h*lu']
cands = lll_rel('swap', w3pure)
for maxc,res,v in cands[:2]:
    print("swap: maxc=%.3e resid=%s"%(mp.mpf(maxc), mp.nstr(res,3)))
