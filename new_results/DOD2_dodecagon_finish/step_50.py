import mpmath as mp, time
mp.mp.dps = 150
t0=time.time()
r3 = mp.sqrt(3)
X = {1: 1-2*r3/3, 2: 1-r3, 5: 1-r3/3, 6: 1-r3/2}
IDX=[1,2,5,6]
Li2 = lambda z: mp.polylog(2,z); Li3 = lambda z: mp.polylog(3,z)
def L1abs(kap,u):
    if kap*u <= 1:
        return -Li2(kap*u)
    else:
        return mp.log(kap)*mp.log(u) + mp.log(u)**2/2 + Li2(1/(kap*u))
def I2abs(kap,u):
    if kap*u <= 1:
        return -mp.log(u)*Li2(kap*u) + Li3(kap*u)
    else:
        return mp.log(kap)*mp.log(u)**2/2 + mp.log(u)**3/3 + mp.log(u)*Li2(1/(kap*u)) + Li3(1/(kap*u))
def Cdef(al,be,u):
    return mp.quad(lambda s: mp.log(1-al*s)*mp.log(1-be*s)/s, [0,u])
def Mc(c, xi):
    kap = c/(c-1)
    lo, hi = 1-xi, mp.mpf(1)
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
for (i,j,k) in [(2,5,6),(1,5,6),(2,1,6)]:
    qd = mp.quad(lambda tt: mp.log(1-X[j]*tt)*mp.log(1-X[k]*tt)/(1-X[i]*tt), [0,1])
    cf = T3a_cf(i,j,k)
    print((i,j,k), "diff:", mp.nstr(abs(qd-cf),4), "cf=", mp.nstr(cf,8), "qd=", mp.nstr(qd,8))
print("elapsed", time.time()-t0)
