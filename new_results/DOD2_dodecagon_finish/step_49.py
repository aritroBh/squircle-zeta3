code = r'''
import mpmath as mp, pickle, time
import sympy as sp
mp.mp.dps = 620
t0=time.time()
D1 = '../DOD1_dodecagon/'
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
X = atoms['X']; T3a = atoms['T3a']
g = pickle.load(open('gdswaps620.pkl','rb'))
SwG, SwD = g['SwG'], g['SwD']
Li2R = lambda z: mp.re(mp.polylog(2, mp.mpc(z)))
Li3R = lambda z: mp.re(mp.polylog(3, mp.mpc(z)))
def L1abs(kap,u):
    if kap*u <= 1:
        return -Li2R(kap*u)
    else:
        return mp.log(kap)*mp.log(u) + mp.log(u)**2/2 + Li2R(1/(kap*u))
def I2abs(kap,u):
    if kap*u <= 1:
        return -mp.log(u)*Li2R(kap*u) + Li3R(kap*u)
    else:
        return mp.log(kap)*mp.log(u)**2/2 + mp.log(u)**3/3 + mp.log(u)*Li2R(1/(kap*u)) + Li3R(1/(kap*u))
def Cdef(al,be,u):
    return mp.quad(lambda s: mp.log(1-al*s)*mp.log(1-be*s)/s, [0,u])
def Mc(c, xi):
    kap = c/(c-1); lo, hi = 1-xi, mp.mpf(1)
    return ( mp.log(abs(1-c))*(mp.log(hi)**2-mp.log(lo)**2)/2
             + I2abs(kap,hi) - I2abs(kap,lo)
             - (mp.log(hi)**3-mp.log(lo)**3)/3 )
def Qab(a,b,xi):
    ka, kb = a/(a-1), b/(b-1)
    A, B = mp.log(abs(1-a)), mp.log(abs(1-b))
    lo, hi = 1-xi, mp.mpf(1)
    def F(u):
        return ( A*mp.log(u) + A*L1abs(kb,u) - A*mp.log(u)**2/2
               + B*L1abs(ka,u) - B*mp.log(u)**2/2
               + Cdef(ka,kb,u)
               - I2abs(ka,u) - I2abs(kb,u)
               + mp.log(u)**3/3 )
    return F(hi)-F(lo)
def T3a_cf(i,j,k):
    xi = X[i]; a = X[j]/xi; b = X[k]/xi
    J = -mp.log(1-xi)**3/3 + Mc(a,xi) + Mc(b,xi) + Qab(a,b,xi)
    return J/xi
i,j,k = 2,5,6
cf_ = T3a_cf(i,j,k)
qd = T3a[(i,j,k)]
print("|T3a_cf - T3a_quad| =", mp.nstr(abs(cf_-qd),4), "t=%ds"%(time.time()-t0))
R_cf = cf_ - (SwG[(i,j,k)]-SwD[(i,j,k)])/X[i]
R_qd = qd - (SwG[(i,j,k)]-SwD[(i,j,k)])/X[i]
print("R_cf =", mp.nstr(R_cf,15))
print("R_qd =", mp.nstr(R_qd,15))
print("|R_cf - R_qd| =", mp.nstr(abs(R_cf-R_qd),4))
# test R_cf against basis via LLL
B = pickle.load(open('basis102_620.pkl','rb'))
names, vals = B['names'], B['vals']
ALL = vals + [R_cf]; AN = names + ['Rcf']
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
print("Rcf-involving vectors (sorted by resid):")
for resid,nrm,vc in hits[:5]:
    print("   maxc=%s resid=%s truefloor=%s" % (mp.nstr(mp.mpf(nrm),3), mp.nstr(resid,3), mp.nstr(mp.mpf(nrm)*mp.mpf(10)**-615,3)))
    if resid < mp.mpf(nrm)*mp.mpf(10)**-600:
        print("      ", [(AN[q],vc[q]) for q in range(N) if vc[q]!=0][:30])
print("elapsed", time.time()-t0)
'''
open('step_13.py','w').write(code)
exec(code)
