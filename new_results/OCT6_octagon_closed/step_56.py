import mpmath as mp, pickle
from flint import fmpz_mat
from fractions import Fraction as Fr
mp.mp.dps = 745
r2=mp.sqrt(2); k1=r2-1; k2=r2/2-1; x1=-k1; x2=-k2
l2=mp.log(2); lu=mp.log(1+r2); z2=mp.pi**2/6; z3=mp.zeta(3)
def AA_of(xv):
    N = int(750/mp.log10(1/abs(xv)))+50
    s=mp.mpf(0); Hn=mp.mpf(0); xp=mp.mpf(1)
    for nn in range(1,N+1):
        Hn+=1/mp.mpf(nn); xp*=xv; s+=xp*Hn/nn**2
    return s
vals = pickle.load(open('basis745.pkl','rb'))
vals['AA_r2h'] = AA_of(r2/2)
vals['AA_sum'] = vals['AA_1']+vals['AA_2']+vals['AA_r2h']
print("A(x1)+A(x2)+A(r2h) =", mp.nstr(vals['AA_sum'],30))
# LLL against w3 basis
w3 = ['z3','Li3_mr2h','Li3_x1','Li3_x2','Li3_k1','Li3_d','Li3_r2h','z2l2','z2lu','l2^3','l2^2lu','l2lu^2','lu^3',
      'Li2_mr2h*l2','Li2_mr2h*lu','Li2_x1*l2','Li2_x1*lu','Li2_x2*l2','Li2_x2*lu','Li2_k1*l2','Li2_k1*lu',
      'Li2_d*l2','Li2_d*lu','Li2_r2h*l2','Li2_r2h*lu']
nb=len(w3); n=2*nb+1; C=mp.mpf(10)**600
vec=[vals['AA_sum']]
for b in w3: vec += [vals[b], r2*vals[b]]
rows=[]
for i in range(n):
    row=[0]*n+[int(mp.nint(C*vec[i]))]; row[i]=1; rows.append(row)
M=fmpz_mat(rows).lll()
cands=[]
for i in range(n):
    v=[int(M[i,j]) for j in range(n)]
    if v[0]==0: continue
    res=abs(sum(mp.mpf(v[j])*vec[j] for j in range(n))); maxc=max(abs(xx) for xx in v)
    cands.append((maxc,res,v))
cands.sort(key=lambda z:z[0])
for maxc,res,v in cands[:2]:
    print("maxc=%.2e resid=%s"%(mp.mpf(maxc),mp.nstr(res,3)))
    if res<mp.mpf('1e-400'):
        for i,b in enumerate(w3):
            a,bb=v[1+2*i],v[2+2*i]
            if a or bb: print("    ",Fr(a,-v[0]),"+",Fr(bb,-v[0]),"r2 *",b)
