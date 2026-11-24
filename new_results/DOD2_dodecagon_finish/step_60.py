import mpmath as mp, pickle
mp.mp.dps = 620
D1 = '../DOD1_dodecagon/'
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
i,j,k = 2,5,6
xi = X[i]; a = X[j]/xi; b = X[k]/xi
ka, kb = a/(a-1), b/(b-1)
A, B = mp.log(abs(1-a)), mp.log(abs(1-b))
lo, hi = 1-xi, mp.mpf(1)
def Mcf(c):
    kap = c/(c-1)
    return ( mp.log(abs(1-c))*(mp.log(hi)**2-mp.log(lo)**2)/2 + I2abs(kap,hi) - I2abs(kap,lo)
             - (mp.log(hi)**3-mp.log(lo)**3)/3 )
# QnonC = (F(1)-F(lo)) - Cre(ka,kb,1) + Cre(ka,kb,lo)
def FnonC(u):
    return ( A*B*mp.log(u) + A*L1abs(kb,u) - A*mp.log(u)**2/2
           + B*L1abs(ka,u) - B*mp.log(u)**2/2
           - I2abs(ka,u) - I2abs(kb,u) + mp.log(u)**3/3 )
QnonC = FnonC(hi)-FnonC(lo)
# swap-remainder pieces
da, db = ka*lo, kb*lo
SwGrem = -mp.log(abs(1-ka))*(Li2R(kb)-Li2R(kb/ka)) + Li3R(kb)
SwDrem = -mp.log(abs(1-da))*(Li2R(db)-Li2R(db/da)) + Li3R(db)
bnd = -mp.log(1-xi)**3/3
R_explicit = (bnd + Mcf(a) + Mcf(b) + QnonC + SwGrem - SwDrem)/xi
R_quad = T3a[(i,j,k)] - (SwG[(i,j,k)]-SwD[(i,j,k)])/xi
print("R_explicit =", mp.nstr(R_explicit, 20))
print("R_quad     =", mp.nstr(R_quad, 20))
print("|diff| =", mp.nstr(abs(R_explicit-R_quad),4))
# Now PSLQ R_explicit against B54_ind quickly at this precision? First try small pslq
import pickle
Q = pickle.load(open('basis_quotient.pkl','rb'))
Bp = pickle.load(open('basis102_620.pkl','rb'))
keep54 = [i for i in Q['keep'] if i < 54]
Bn = [Bp['names'][q] for q in keep54]; Bv = [Bp['vals'][q] for q in keep54]
rel = mp.pslq([R_quad]+Bv, tol=mp.mpf('1e-560'), maxcoeff=10**7, maxsteps=10000)
print("pslq:", "None" if rel is None else rel)
if rel:
    resid = rel[0]*R_quad + sum(rel[q+1]*Bv[q] for q in range(len(Bv)))
    print("resid:", mp.nstr(abs(resid),4))
    print([(Bn[q],rel[q+1]) for q in range(len(Bv)) if rel[q+1]!=0])
