code = r'''
import mpmath as mp, pickle, time
import sympy as sp
mp.mp.dps = 620
t0 = time.time()
D1 = '../DOD1_dodecagon/'
r3m = mp.sqrt(3)
r3s = sp.sqrt(3)
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
X = atoms['X']; Lam = atoms['Lam']; T3a = atoms['T3a']
g = pickle.load(open('gdswaps620.pkl','rb'))
SwG, SwD = g['SwG'], g['SwD']
cf = pickle.load(open('atom_coefs_pq.pkl','rb'))
B = pickle.load(open('basis102_620.pkl','rb'))
Q = pickle.load(open('basis_quotient.pkl','rb'))
names, vals = B['names'], B['vals']
keep = Q['keep']
Bind_names = [names[i] for i in keep]; Bind_vals = [vals[i] for i in keep]
# 1/x_i in Q(sqrt3)
invX = {1: (-3,-2), 2: (-sp.Rational(1,2),-sp.Rational(1,2)), 5: (sp.Rational(3,2),sp.Rational(1,2)), 6: (4,2)}
for i,(u,v) in invX.items():
    assert abs((u+v*r3m)*X[i]-1) < mp.mpf('1e-600'), i
# Build Xi combo coefficients
IDX=[1,2,5,6]
Xi_p = mp.mpf(0); Xi_q = mp.mpf(0); Xi_full = mp.mpf(0)
terms = []  # (kind, key, p, q)
for (fam,ij) in cf:
    if fam=='Lam':
        p,q = cf[(fam,ij)]
        terms.append(('Lam', ij, p, q))
    elif fam=='T3a':
        P,Qq = cf[(fam,ij)]
        i = ij[0]
        u,v = invX[i]
        Pp = P*u + 3*Qq*v
        Qp = P*v + Qq*u
        terms.append(('SwGD', ij, Pp, Qp))
Xi_terms = []
for (kind,key,p,q) in terms:
    if kind=='Lam':
        v = Lam[key]
    else:
        v = SwG[key]-SwD[key]
    p = sp.Rational(p); q = sp.Rational(q)
    Xi_p += mp.mpf(str(p))*v
    Xi_q += mp.mpf(str(q))*v
    Xi_full += (mp.mpf(str(p))+mp.mpf(str(q))*r3m)*v
    Xi_terms.append((kind,key,str(p),str(q)))
print("Xi_p =", mp.nstr(Xi_p, 20))
print("Xi_q =", mp.nstr(Xi_q, 20))
print("|Xi_full - (Xi_p + sqrt3 Xi_q)| =", mp.nstr(abs(Xi_full-(Xi_p+r3m*Xi_q)),4))
# direct remainder from S: Xi_direct = sum cLam*Lam + sum cT3a*T3a
XiD_p = mp.mpf(0); XiD_q = mp.mpf(0)
for (fam,ij) in cf:
    if fam=='Lam':
        p,q = cf[(fam,ij)]; v = Lam[ij]
    elif fam=='T3a':
        p,q = cf[(fam,ij)]; v = T3a[ij]
    else:
        continue
    XiD_p += mp.mpf(str(sp.Rational(p)))*v
    XiD_q += mp.mpf(str(sp.Rational(q)))*v
print("XiD_p =", mp.nstr(XiD_p, 20))
print("XiD_q =", mp.nstr(XiD_q, 20))
print("|Xi_p - XiD_p| =", mp.nstr(abs(Xi_p-XiD_p),4), " |Xi_q - XiD_q| =", mp.nstr(abs(Xi_q-XiD_q),4))
pickle.dump({'Xi_p':Xi_p,'Xi_q':Xi_q,'XiD_p':XiD_p,'XiD_q':XiD_q,'Xi_terms':Xi_terms},
            open('Xi_620.pkl','wb'))
# targets: 12 R's + Xi_p, Xi_q, XiD_p, XiD_q
tvals=[]; tnames=[]
for (i,j,k) in sorted(T3a.keys()):
    tvals.append(T3a[(i,j,k)] - (SwG[(i,j,k)]-SwD[(i,j,k)])/X[i]); tnames.append('R%d_%d%d'%(i,j,k))
tvals += [Xi_p, Xi_q, XiD_p, XiD_q]; tnames += ['Xi_p','Xi_q','XiD_p','XiD_q']
ALL = Bind_vals + tvals
AN = Bind_names + tnames
N = len(ALL)
print("lattice size:", N, "t=%ds"%(time.time()-t0))
C = 10**540
rows = []
for i in range(N):
    rr_ = [0]*N + [int(mp.nint(C*ALL[i]))]
    rr_[i] = 1
    rows.append(rr_)
M = sp.Matrix(rows)
print("LLL start t=%ds"%(time.time()-t0))
BM = M.lll(delta=sp.Rational(3,4))
print("LLL done t=%ds"%(time.time()-t0))
NB = len(Bind_vals)
found = {}; bas = 0
for i in range(BM.rows):
    vc = [int(BM[i,j]) for j in range(N)]
    nrm = max(abs(x) for x in vc)
    if nrm==0 or nrm > 10**6: continue
    resid = abs(sum(mp.mpf(vc[q])*ALL[q] for q in range(N)))
    if resid < mp.mpf('1e-500'):
        sup_t = [AN[q] for q in range(NB,N) if vc[q]!=0]
        if sup_t:
            for nm in sup_t: found.setdefault(nm,[]).append((vc,resid,nrm))
        else:
            bas += 1
print("leftover pure-basis rels:", bas)
for nm in tnames:
    if nm in found:
        vc,resid,nrm = sorted(found[nm], key=lambda t:t[2])[0]
        sup = [(AN[q],vc[q]) for q in range(N) if vc[q]!=0]
        print("REDUCES:", nm, "maxc=%d resid=%s"%(nrm,mp.nstr(resid,3)))
        print("   ", sup[:20], ("+%d"%(len(sup)-20) if len(sup)>20 else ""))
    else:
        print("NO REL:", nm)
pickle.dump({'found':{k:[(vc,str(r),n) for vc,r,n in v] for k,v in found.items()}}, open('xi_test620.pkl','wb'))
print("elapsed", time.time()-t0)
'''
open('step_07.py','w').write(code)
exec(code)
