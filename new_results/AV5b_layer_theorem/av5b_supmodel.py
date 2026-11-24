
import mpmath as mp
mp.mp.dps = 55
import json, sys

def geom(m):
    rho = mp.tan(mp.pi/(2*m))
    v0 = (mp.mpf(0), rho)
    v1 = (rho*mp.sin(mp.pi/m), rho*mp.cos(mp.pi/m))
    th0 = mp.pi/2 - mp.pi/(2*m)
    N = (mp.cos(th0), mp.sin(th0))
    H = mp.sin(mp.pi/(2*m))
    dv = (v1[0]-v0[0], v1[1]-v0[1])
    L = mp.sqrt(dv[0]**2+dv[1]**2)
    ns = [(mp.cos(i*mp.pi/m), mp.sin(i*mp.pi/m)) for i in range(1, m+1)]
    return v0, v1, th0, N, H, L, dv, ns

class Z:
    def __init__(self, m, full):
        self.m, self.full = m, full
        self.v0, self.v1, self.th0, self.N, self.H, self.L, self.dv, self.ns = geom(m)
        self.j = {m-1, m-2} if not full else set(range(m))
    def ab(self, t):
        v0, dv, ns, N = self.v0, self.dv, self.ns, self.N
        xt = (v0[0]+t*dv[0], v0[1]+t*dv[1])
        a = [n[0]*xt[0]+n[1]*xt[1] for n in ns]
        b = [n[0]*N[0]+n[1]*N[1] for n in ns]
        ap = [n[0]*dv[0]+n[1]*dv[1] for n in ns]
        return a, b, ap
    def F(self, t, z, p):
        a, b, ap = self.ab(t)
        s = mp.mpf(0)
        for i in range(self.m):
            if self.full or i in self.j:
                aa = a[i]+z*b[i]
                s += abs(aa)**p
            else:
                s += abs(a[i])
        return s - 1
    def Fz(self, t, z, p):
        a, b, ap = self.ab(t)
        s = mp.mpf(0)
        for i in range(self.m):
            if self.full or i in self.j:
                aa = a[i]+z*b[i]
                if aa != 0: s += p*b[i]*mp.sign(aa)*abs(aa)**(p-1)
        return s
    def Ft(self, t, z, p):
        a, b, ap = self.ab(t)
        s = mp.mpf(0)
        for i in range(self.m):
            if self.full or i in self.j:
                aa = a[i]+z*b[i]
                if aa != 0: s += p*mp.sign(aa)*abs(aa)**(p-1)*ap[i]
            else:
                s += mp.sign(a[i])*ap[i]
        return s
    def zs(self, t, p):
        if t <= 0 or t >= 1: return mp.mpf(0)
        a, b, ap = self.ab(t)
        S1 = mp.fsum(abs(ai)*mp.log(abs(ai)) if (self.full or i in self.j) and ai != 0 else mp.mpf(0)
                     for i, ai in enumerate(a))
        Sb = mp.fsum(mp.sign(a[i])*b[i] for i in range(self.m) if (self.full or i in self.j))
        z = -(p-1)*S1/Sb if Sb != 0 else mp.mpf(0)
        for _ in range(80):
            fz = self.F(t, z, p)
            if abs(fz) < mp.mpf('1e-48'): break
            dz = fz/self.Fz(t, z, p)
            z -= dz
            if abs(dz) < mp.mpf('1e-48'): break
        return z
    def zp(self, t, p):
        z = self.zs(t, p)
        return -self.Ft(t, z, p)/self.Fz(t, z, p)

def h_loc(Zo, u, p):
    # support of model curve at angle th0 + u
    tgt = -Zo.L*mp.tan(u)
    f = lambda t: Zo.zp(t, p) - tgt
    if abs(u) < mp.mpf('1e-45'):
        t = mp.mpf('0.5')
    else:
        lo, hi = mp.mpf('1e-28'), 1-mp.mpf('1e-28')
        flo, fhi = f(lo), f(hi)
        for _ in range(300):
            mid = (lo+hi)/2
            fm = f(mid)
            if abs(fm) < mp.mpf('1e-42'): break
            if (fm > 0) == (flo > 0):
                lo, flo = mid, fm
            else:
                hi, fhi = mid, fm
        t = (lo+hi)/2
    z = Zo.zs(t, p)
    xt0 = Zo.v0[0]+t*Zo.dv[0]; xt1 = Zo.v0[1]+t*Zo.dv[1]
    mx = (Zo.v0[0]+Zo.v1[0])/2; my = (Zo.v0[1]+Zo.v1[1])/2
    Tdir = (Zo.dv[0]/Zo.L, Zo.dv[1]/Zo.L)
    ell = Tdir[0]*(xt0-mx)+Tdir[1]*(xt1-my)
    return (Zo.H+z)*mp.cos(u) + ell*mp.sin(u)

def dP(Zo, p):
    W = mp.pi/(2*Zo.m)
    f = lambda u: h_loc(Zo, u, p) - (Zo.H*mp.cos(u) + (Zo.L/2)*abs(mp.sin(u)))
    return 2*Zo.m*mp.quad(f, [-W, 0, W])

if __name__ == "__main__":
    mm = int(sys.argv[1]); full = int(sys.argv[2])
    Zo = Z(mm, full)
    out = []
    for dstr in sys.argv[3].split(","):
        p = 1+mp.mpf(dstr)
        v = dP(Zo, p)
        out.append((dstr, v))
        print(f"m={mm} full={full} d={dstr}: dP={mp.nstr(v,35)}", flush=True)
    json.dump([[d, str(v)] for d, v in out], open(f"ckpt_sup_m{mm}_f{full}.json","w"))
