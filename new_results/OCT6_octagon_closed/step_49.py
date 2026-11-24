import mpmath as mp, pickle
from flint import fmpz_mat
from fractions import Fraction as Fr
mp.mp.dps = 745
vals = pickle.load(open('basis745.pkl','rb'))
r2m=mp.sqrt(2)
vals['Li3b_oct4'] = mp.polylog(3, (2+r2m)/4)
vals['Li3c_oct4'] = mp.polylog(3, 12*r2m-16)
w3 = ['z3','Li3_mr2h','Li3_x1','Li3_x2','Li3_k1','Li3_d','Li3_r2h','z2l2','z2lu','l2^3','l2^2lu','l2lu^2','lu^3']
def find(tgt):
    nb=len(w3); n=2*nb+1; C=mp.mpf(10)**600
    vec=[vals[tgt]]
    for b in w3: vec += [vals[b], r2m*vals[b]]
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
    cands.sort(key=lambda z:z[0])
    maxc,res,v = cands[0]
    out=[(Fr(v[1+2*i],-v[0]),Fr(v[2+2*i],-v[0]),b) for i,b in enumerate(w3) if v[1+2*i] or v[2+2*i]]
    return maxc,res,out
for nm in ['Li3b_oct4','Li3c_oct4']:
    maxc,res,out = find(nm)
    print(nm, "maxc",maxc,"resid", mp.nstr(res,3))
    for a,b,bb in out: print("   ", a, "+", b, "r2 *", bb)
