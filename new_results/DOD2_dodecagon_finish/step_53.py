import mpmath as mp, time
mp.mp.dps = 150
t0=time.time()
r3 = mp.sqrt(3)
X = {1: 1-2*r3/3, 2: 1-r3, 5: 1-r3/3, 6: 1-r3/2}
IDX=[1,2,5,6]
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
def Cre(al,be,u):
    # int_0^u ln|1-al s| ln|1-be s|/s ds, splitting at interior singularities
    pts = [mp.mpf(0)]
    for kap in (al,be):
        if kap > 0:
            s0 = 1/kap
            if 0 < s0 < u: pts.append(s0)
    pts.append(u)
    pts = sorted(set(pts))
    tot = mp.mpf(0)
    for a_,b_ in zip(pts[:-1], pts[1:]):
        tot += mp.quad(lambda s: mp.log(abs(1-al*s))*mp.log(abs(1-be*s))/s, [a_,b_])
    return tot
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
        return ( A*B*mp.log(u) + A*L1abs(kb,u) - A*mp.log(u)**2/2
               + B*L1abs(ka,u) - B*mp.log(u)**2/2
               + Cre(ka,kb,u)
               - I2abs(ka,u) - I2abs(kb,u)
               + mp.log(u)**3/3 )
    return F(hi)-F(lo)
def T3a_cf(i,j,k):
    xi = X[i]; a = X[j]/xi; b = X[k]/xi
    J = -mp.log(1-xi)**3/3 + Mc(a,xi) + Mc(b,xi) + Qab(a,b,xi)
    return J/xi
worst = mp.mpf(0)
for i in IDX:
    for j in IDX:
        for k in IDX:
            if len({i,j,k})==3 and j<k:
                qd = mp.quad(lambda tt: mp.log(1-X[j]*tt)*mp.log(1-X[k]*tt)/(1-X[i]*tt), [0,1])
                cf = T3a_cf(i,j,k)
                d = abs(qd-cf); worst = max(worst,d)
                print((i,j,k), "diff:", mp.nstr(d,4))
print("worst:", mp.nstr(worst,4), "elapsed", time.time()-t0)
