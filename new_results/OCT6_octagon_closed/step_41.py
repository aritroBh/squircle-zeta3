import pickle, mpmath as mp
from flint import fmpz_mat
mp.mp.dps = 745
vals = pickle.load(open('basis745.pkl','rb'))
r2m = mp.sqrt(2)
# add products needed
vals['l2^2']=vals['l2']**2; vals['l2*lu']=vals['l2']*vals['lu']; vals['lu^2']=vals['lu']**2
named_w2 = ['z2','l2^2','l2*lu','lu^2','Li2_x1','Li2_x2','Li2_k1','Li2_d','Li2_r2h','Li2_mr2h']
named_w3 = ['z3','Li3_x1','Li3_x2','Li3_k1','Li3_d','Li3_r2h','Li3_mr2h','z2l2','z2lu','l2^3','l2^2lu','l2lu^2','lu^3',
 'Li2_x1*l2','Li2_x2*l2','Li2_k1*l2','Li2_d*l2','Li2_r2h*l2','Li2_mr2h*l2',
 'Li2_x1*lu','Li2_x2*lu','Li2_k1*lu','Li2_d*lu','Li2_r2h*lu','Li2_mr2h*lu']
def lll_rel(target, basis, dps_scale=0.80):
    nb=len(basis); n=2*nb+1
    C=mp.mpf(10)**int(mp.mp.dps*dps_scale)
    vec=[vals[target]]
    for b in basis: vec += [vals[b], r2m*vals[b]]
    rows=[]
    for i in range(n):
        row=[0]*n+[int(mp.nint(C*vec[i]))]; row[i]=1; rows.append(row)
    M=fmpz_mat(rows).lll()
    best=None
    for i in range(n):
        v=[int(M[i,j]) for j in range(n)]
        if v[0]==0: continue
        res=abs(sum(mp.mpf(v[j])*vec[j] for j in range(n))); maxc=max(abs(x) for x in v)
        if best is None or maxc<best[0]: best=(maxc,res,v)
    return best
def greedy(named):
    basis=[]; rels=[]
    for c in named:
        if not basis: basis.append(c); continue
        maxc,res,v = lll_rel(c, basis)
        if res < mp.mpf('1e-600') and maxc < 1e6:
            from fractions import Fraction as Fr
            rel = [(basis[i], Fr(v[1+2*i],-v[0]), Fr(v[2+2*i],-v[0])) for i in range(len(basis)) if v[1+2*i] or v[2+2*i]]
            rels.append((c, rel, res, maxc))
        else:
            basis.append(c)
    return basis, rels
b2, r2rels = greedy(named_w2)
print("W2 independent basis (%d of %d):"%(len(b2),len(named_w2)), b2)
for c,rel,res,maxc in r2rels:
    print(f"  REL {c} = ", [(nm,str(a),str(b)) for nm,a,b in rel], "resid", mp.nstr(res,3), "maxc", maxc)
b3, r3rels = greedy(named_w3)
print("\nW3 independent basis (%d of %d):"%(len(b3),len(named_w3)), b3)
for c,rel,res,maxc in r3rels:
    print(f"  REL {c} = ", [(nm,str(a),str(b)) for nm,a,b in rel], "resid", mp.nstr(res,3), "maxc", maxc)
pickle.dump({'b2':b2,'r2rels':r2rels,'b3':b3,'r3rels':r3rels}, open('indep.pkl','wb'))
