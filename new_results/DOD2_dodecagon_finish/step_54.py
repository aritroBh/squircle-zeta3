code = r'''
import mpmath as mp, pickle, time
import sympy as sp
mp.mp.dps = 620
t0=time.time()
D1 = '../DOD1_dodecagon/'
r3 = mp.sqrt(3)
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
X = atoms['X']; T3a = atoms['T3a']
g = pickle.load(open('gdswaps620.pkl','rb'))
SwG, SwD = g['SwG'], g['SwD']
Li2R = lambda z: mp.re(mp.polylog(2, mp.mpc(z)))
Li3R = lambda z: mp.re(mp.polylog(3, mp.mpc(z)))
def L1abs(kap,u):
    if kap*u <= 1: return -Li2R(kap*u)
    else: return mp.log(kap)*mp.log(u) + mp.log(u)**2/2 + Li2R(1/(kap*u))
def I2abs(kap,u):
    if kap*u <= 1: return -mp.log(u)*Li2R(kap*u) + Li3R(kap*u)
    else: return mp.log(kap)*mp.log(u)**2/2 + mp.log(u)**3/3 + mp.log(u)*Li2R(1/(kap*u)) + Li3R(1/(kap*u))
def Cre(al,be,u):
    pts = [mp.mpf(0)]
    for kap in (al,be):
        if kap > 0:
            s0 = 1/kap
            if 0 < s0 < u: pts.append(s0)
    pts.append(u); pts = sorted(set(pts))
    tot = mp.mpf(0)
    for a_,b_ in zip(pts[:-1], pts[1:]):
        tot += mp.quad(lambda s: mp.log(abs(1-al*s))*mp.log(abs(1-be*s))/s, [a_,b_])
    return tot
def Mc(c, xi):
    kap = c/(c-1); lo, hi = 1-xi, mp.mpf(1)
    return ( mp.log(abs(1-c))*(mp.log(hi)**2-mp.log(lo)**2)/2 + I2abs(kap,hi) - I2abs(kap,lo)
             - (mp.log(hi)**3-mp.log(lo)**3)/3 )
def Qab(a,b,xi):
    ka, kb = a/(a-1), b/(b-1)
    A, B = mp.log(abs(1-a)), mp.log(abs(1-b))
    lo, hi = 1-xi, mp.mpf(1)
    def F(u):
        return ( A*B*mp.log(u) + A*L1abs(kb,u) - A*mp.log(u)**2/2
               + B*L1abs(ka,u) - B*mp.log(u)**2/2 + Cre(ka,kb,u)
               - I2abs(ka,u) - I2abs(kb,u) + mp.log(u)**3/3 )
    return F(hi)-F(lo)
def T3a_cf(i,j,k):
    xi = X[i]; a = X[j]/xi; b = X[k]/xi
    J = -mp.log(1-xi)**3/3 + Mc(a,xi) + Mc(b,xi) + Qab(a,b,xi)
    return J/xi
i,j,k = 2,5,6
cf_ = T3a_cf(i,j,k)
qd = T3a[(i,j,k)]
print("|T3a_cf - quad| =", mp.nstr(abs(cf_-qd),4), "t=%ds"%(time.time()-t0))
R = qd - (SwG[(i,j,k)]-SwD[(i,j,k)])/X[i]
print("R(2,5,6) =", mp.nstr(R,15))
B = pickle.load(open('basis102_620.pkl','rb'))
names, vals = B['names'], B['vals']
ALL = vals[:54] + [R]; AN = names[:54] + ['R']
N = len(ALL)
C = 10**540
rows = []
for q in range(N):
    rr_ = [0]*N + [int(mp.nint(C*ALL[q]))]; rr_[q] = 1
    rows.append(rr_)
BM = sp.Matrix(rows).lll(delta=sp.Rational(3,4))
hits=[]
for r in range(BM.rows):
    vc = [int(BM[r,q]) for q in range(N)]
    if vc[N-1]==0: continue
    nrm = max(abs(x) for x in vc)
    resid = abs(sum(mp.mpf(vc[q])*ALL[q] for q in range(N)))
    hits.append((resid, nrm, vc))
hits.sort()
for resid,nrm,vc in hits[:4]:
    print("   maxc=%s resid=%s truefloor=%s" % (mp.nstr(mp.mpf(nrm),3), mp.nstr(resid,3), mp.nstr(mp.mpf(nrm)*mp.mpf(10)**-615,3)))
    if resid < mp.mpf(nrm)*mp.mpf(10)**-600:
        print("      ", [(AN[q],vc[q]) for q in range(N) if vc[q]!=0][:26])
print("elapsed", time.time()-t0)
'''
open('step_14.py','w').write(code)
exec(code)
