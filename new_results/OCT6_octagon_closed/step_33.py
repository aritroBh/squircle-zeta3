import pickle, mpmath as mp
from flint import fmpz_mat
from fractions import Fraction as Fr
mp.mp.dps = 745
vals = pickle.load(open('basis745.pkl','rb'))
r2 = mp.sqrt(2)
def lll_rel(target, basis_names, dps_scale=0.75):
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
def show(v, basis_names):
    out=[]
    for i,b in enumerate(basis_names):
        a,bb=v[1+2*i],v[2+2*i]
        if a or bb:
            ca, cb = Fr(a,-v[0]), Fr(bb,-v[0])
            out.append((ca,cb,b))
    return out
w3pure = ['z3','Li3_mr2h','Li3_x1','Li3_x2','Li3_k1','Li3_d','Li3_r2h','z2l2','z2lu','l2^3','l2^2lu','l2lu^2','lu^3',
      'Li2_mr2h*l2','Li2_mr2h*lu','Li2_x1*l2','Li2_x1*lu','Li2_x2*l2','Li2_x2*lu','Li2_k1*l2','Li2_k1*lu',
      'Li2_d*l2','Li2_d*lu','Li2_r2h*l2','Li2_r2h*lu']
for tgt in ['Lam','swap','AA_2']:
    cands = lll_rel(tgt, w3pure)
    print("="*25, tgt)
    for maxc,res,v in cands[:2]:
        print("  maxc=%.3e resid=%s"%(mp.mpf(maxc), mp.nstr(res,3)))
        if res < mp.mpf('1e-400'):
            for ca,cb,b in show(v, w3pure): print("    %+s %+s*r2  * %s"%(ca,cb,b))
