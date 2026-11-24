
"""PG5 step_05: letter group of the 14-gon. Letters {2, a=2cos(2pi/7), u=2cos(pi/7)=a^2+a-1}.
factor_ft: x -> (sign, e2, ea, eu) with x = +- 2^e2 a^ea u^eu, via norm + embedding lattice.
Verifies exactly by multiplying back. Then factor all geometry numbers, all pairwise
differences, all gamma/delta cross points, and the full CR set of {0,1,inf,p1,p2,p5,p6}.
"""
import mpmath as mp, pickle, sys
sys.path.insert(0,'.')
import ftfield as F
mp.mp.dps = 80
A = F.AL
U1 = F.ft(-1,1,1)           # u = a^2+a-1 = 2cos(pi/7)
TWO = F.ft(2)
EMB = F.get_EMB()
def emb(v,k):
    A0,B,C,d = v
    a = EMB[k]
    return (mp.mpf(int(A0.numerator))/A0.denominator + mp.mpf(int(B.numerator))/B.denominator*a + mp.mpf(int(C.numerator))/C.denominator*a*a)/d
la = [mp.log(abs(emb(A,k))) for k in range(3)]
lu = [mp.log(abs(emb(U1,k))) for k in range(3)]
def factor_ft(x):
    x = F.norm_t(x)
    if x[:3]==(F.Fr(0),)*3: return None
    sgn = 1
    if emb(x,0) < 0: sgn = -1
    N = F.fnorm(x)
    # only prime 2 can appear (inert); check
    num, den = abs(N.numerator), N.denominator
    e2n = 0
    while num%2==0 and num>0: num//=2; e2n+=1
    e2d = 0
    while den%2==0 and den>0: den//=2; e2d+=1
    if num!=1 or den!=1: return ('NOT_SUNIT', N)
    assert (e2n-e2d)%3==0, ('bad 2-adic valuation', N)
    e2 = (e2n-e2d)//3
    xu = F.fmul(x, F.fpow(TWO, -e2))
    lx = [mp.log(abs(emb(xu,k))) for k in range(3)]
    # solve [la;lu]^T (ea,eu) = lx using rows 0,1, verify row 2
    import numpy as np
    M = np.array([[float(la[0]),float(lu[0])],[float(la[1]),float(lu[1])]])
    b = np.array([float(lx[0]),float(lx[1])])
    sol = np.linalg.solve(M,b)
    ea, eu = int(round(sol[0])), int(round(sol[1]))
    for (da,du) in [(0,0),(1,0),(-1,0),(0,1),(0,-1)]:
        cand = F.fmul(F.fpow(A,ea+da), F.fpow(U1,eu+du))
        if F.feq(cand, xu) or F.feq(cand, F.fneg(xu)):
            return (sgn, e2, ea+da, eu+du)
    return ('UNIT_FAIL', sol)
# --- geometry ---
g = pickle.load(open('geom14.pkl','rb'))
X = g['X']; B = g['B']; A0 = g['A0']; hh = g['h']
names = {'A0':A0,'h':hh,'B1':B[1],'B2':B[2],'B5':B[5],'B6':B[6]}
for i in [1,2,5,6]:
    names[f'x{i}']=X[i]; names[f'1-x{i}']=F.fsub(F.ONE,X[i])
# pairwise differences of x's
for i in [1,2,5,6]:
    for j in [1,2,5,6]:
        if i<j: names[f'x{i}-x{j}']=F.fsub(X[i],X[j])
# gamma/delta points: gam_ij = x_j/(x_j-x_i), del_ij = gam_ij(1-x_i)
gam = {}; delta = {}
for i in [1,2,5,6]:
    for j in [1,2,5,6]:
        if i==j: continue
        gij = F.fdiv(X[j], F.fsub(X[j],X[i]))
        gam[(i,j)] = gij
        delta[(i,j)] = F.fmul(gij, F.fsub(F.ONE,X[i]))
        names[f'g{i}{j}']=gij; names[f'1-g{i}{j}']=F.fsub(F.ONE,gij)
        names[f'd{i}{j}']=delta[(i,j)]; names[f'1-d{i}{j}']=F.fsub(F.ONE,delta[(i,j)])
# p-points and CR set
p = {i: F.fdiv(F.ONE, X[i]) for i in [1,2,5,6]}
for i in [1,2,5,6]:
    names[f'p{i}']=p[i]
    for j in [1,2,5,6]:
        if i<j: names[f'p{i}-p{j}']=F.fsub(p[i],p[j])
# CR set over {0,1,inf,p1,p2,p5,p6}
INF='inf'
pts = {'0':F.ZERO,'1':F.ONE,INF:None,'p1':p[1],'p2':p[2],'p5':p[5],'p6':p[6]}
def dsub(a,b):
    if a==INF or b==INF: return INF
    return F.fsub(pts[a],pts[b])
def crat(w,x,y,z):
    num1, num2 = dsub(w,y), dsub(x,z)
    den1, den2 = dsub(w,z), dsub(x,y)
    inf = 0; fac = []
    for v in (num1,num2):
        if v==INF: inf+=1
        else: fac.append((v,1))
    for v in (den1,den2):
        if v==INF: inf-=1
        else: fac.append((v,-1))
    assert inf==0
    r = F.ONE
    for v,e in fac: r = F.fmul(r, F.fpow(v,e))
    return F.norm_t(r)
import itertools
CRset = {}
pnames = ['0','1',INF,'p1','p2','p5','p6']
for q4 in itertools.combinations(pnames,4):
    for (w,x,y,z) in [(0,1,2,3),(0,1,3,2),(0,2,1,3),(0,2,3,1),(0,3,1,2),(0,3,2,1)]:
        q = [q4[w],q4[x],q4[y],q4[z]]
        cr = crat(*q)
        key = cr
        if key not in CRset:
            CRset[key] = (cr, F.norm_t(F.fsub(F.ONE,cr)))
print("number of distinct CRs (14-gon edge):", len(CRset))
for nm,(cr,oc) in list(CRset.items())[:5]: pass
# add CRs and 1-CRs to names
for k,(cr,oc) in enumerate(CRset.values()):
    names[f'cr{k}']=cr; names[f'1-cr{k}']=oc
# ---- factor everything ----
fails = {}
facs = {}
for nm,v in names.items():
    f = factor_ft(v)
    if isinstance(f,tuple) and f and f[0]=='NOT_SUNIT' or (isinstance(f,tuple) and len(f)==2 and f[0]=='UNIT_FAIL'):
        fails[nm]=f
    else:
        facs[nm]=f
        # verify exactly
        s,e2,ea,eu = f
        back = F.fmul(F.fpow(TWO,e2), F.fmul(F.fpow(A,ea), F.fpow(U1,eu)))
        ok = F.feq(back, F.fmul(v,F.ft(s))) if s==1 else F.feq(back, F.fneg(v))
        assert ok, (nm, f, v)
print("factorization failures:", fails)
print("total factored:", len(facs), "of", len(names))
import collections
expcount = collections.Counter()
for nm,f in facs.items():
    for e in f[1:]:
        if abs(e)>6: expcount[nm]+=1
print("large exponents:", dict(expcount))
pickle.dump({'facs':facs,'CRset':list(CRset.values()),'gam':gam,'delta':delta,'p':p}, open('letters14.pkl','wb'))
open('progress.log','a').write(f'PG5 step_05: letters {{2,a,u}}; CR set size {len(CRset)}; all {len(facs)} factorizations verified exactly; failures: {fails}\n')
