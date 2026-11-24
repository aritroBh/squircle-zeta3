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
i,j,k = 2,5,6
xi = X[i]; a = X[j]/xi; b = X[k]/xi
ka, kb = a/(a-1), b/(b-1)
A, B = mp.log(abs(1-a)), mp.log(abs(1-b))
lo, hi = 1-xi, mp.mpf(1)
# verify each piece against quadrature
# M(a) = int_0^xi La(s) ln(1-s)/(1-s) ds, La = ln((1-a s)/(1-s))
def Mquad(c):
    return mp.quad(lambda s: mp.log(abs((1-c*s)/(1-s)))*mp.log(1-s)/(1-s), [0,xi])
Ma = Mquad(a); Mb = Mquad(b)
def Mcf(c):
    kap = c/(c-1)
    return ( mp.log(abs(1-c))*(mp.log(hi)**2-mp.log(lo)**2)/2 + I2abs(kap,hi) - I2abs(kap,lo)
             - (mp.log(hi)**3-mp.log(lo)**3)/3 )
print("M(a) quad vs cf:", mp.nstr(abs(Ma-Mcf(a)),4), " M(b):", mp.nstr(abs(Mb-Mcf(b)),4))
# Q = int_0^xi La Lb/(1-s) ds
def Qquad():
    return mp.quad(lambda s: mp.log(abs((1-a*s)/(1-s)))*mp.log(abs((1-b*s)/(1-s)))/(1-s), [0,xi])
Qq = Qquad()
def Qcf():
    def F(u):
        return ( A*B*mp.log(u) + A*L1abs(kb,u) - A*mp.log(u)**2/2
               + B*L1abs(ka,u) - B*mp.log(u)**2/2 + Cre(ka,kb,u)
               - I2abs(ka,u) - I2abs(kb,u) + mp.log(u)**3/3 )
    return F(hi)-F(lo)
print("Q quad vs cf:", mp.nstr(abs(Qq-Qcf()),4))
# boundary term
bnd = -mp.log(1-xi)**3/3
bndq = mp.quad(lambda s: mp.log(1-s)**2/(1-s), [0,xi])
print("bnd:", mp.nstr(abs(bnd-bndq),4))
J = bnd + Ma + Mb + Qq
print("J/xi vs T3a_quad:", mp.nstr(abs(J/xi - T3a[(i,j,k)]),4))
