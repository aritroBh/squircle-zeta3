
import mpmath as mp
mp.mp.dps = 50
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

class Loc:
    def __init__(self, m):
        self.m = m
        self.v0, self.v1, self.th0, self.N, self.H, self.L, self.dv, self.ns = geom(m)
        self.j = {m-1, m-2}
        self.alpha = mp.tan(mp.pi/(2*m))
    def ab(self, t):
        v0, dv, ns, N = self.v0, self.dv, self.ns, self.N
        xt = (v0[0]+t*dv[0], v0[1]+t*dv[1])
        a = [n[0]*xt[0]+n[1]*xt[1] for n in ns]
        b = [n[0]*N[0]+n[1]*N[1] for n in ns]
        ap = [n[0]*dv[0]+n[1]*dv[1] for n in ns]
        return a, b, ap
    def F(self, t, z, p):
        a, b, ap = self.ab(t)
        s = abs(a[self.m-1]+z*b[self.m-1])**p + abs(a[self.m-2]+z*b[self.m-2])**p
        for i in range(self.m):
            if i not in self.j: s += abs(a[i])
        return s - 1
    def Fz(self, t, z, p):
        a, b, ap = self.ab(t)
        s = mp.mpf(0)
        for i in self.j:
            aa = a[i]+z*b[i]
            if aa != 0: s += p*b[i]*mp.sign(aa)*abs(aa)**(p-1)
        return s
    def Ft(self, t, z, p):
        a, b, ap = self.ab(t)
        s = mp.mpf(0)
        for i in self.j:
            aa = a[i]+z*b[i]
            if aa != 0: s += p*mp.sign(aa)*abs(aa)**(p-1)*ap[i]
        for i in range(self.m):
            if i not in self.j: s += mp.sign(a[i])*ap[i]
        return s
    def zs(self, t, p):
        if t <= 0 or t >= 1: return mp.mpf(0)
        a, b, ap = self.ab(t)
        S1 = mp.fsum(abs(a[i])*mp.log(max(abs(a[i]), mp.mpf('1e-300'))) for i in self.j)
        Sb = mp.fsum(mp.sign(a[i])*b[i] for i in self.j)
        z = -(p-1)*S1/Sb
        for _ in range(60):
            fz = self.F(t, z, p)
            if abs(fz) < mp.mpf('1e-44'): break
            dz = fz/self.Fz(t, z, p)
            z -= dz
            if abs(dz) < mp.mpf('1e-44'): break
        return z
    def zp(self, t, p):
        z = self.zs(t, p)
        return -self.Ft(t, z, p)/self.Fz(t, z, p)
    def h(self, u, p):
        # support of the local model curve at angle th0 + u
        if abs(u) < mp.mpf('1e-44'):
            t = mp.mpf('0.5')
        else:
            f = lambda t: self.zp(t, p)*mp.cos(u) + self.L*mp.sin(u)
            lo, hi = mp.mpf('1e-30'), 1-mp.mpf('1e-30')
            flo, fhi = f(lo), f(hi)
            if flo <= 0:
                t = lo
            elif fhi >= 0:
                t = hi
            else:
                for _ in range(150):
                    mid = (lo+hi)/2
                    fm = f(mid)
                    if abs(fm) < mp.mpf('1e-40'): break
                    if (fm>0)==(flo>0): lo, flo = mid, fm
                    else: hi = mid
                t = (lo+hi)/2
        z = self.zs(t, p)
        xt0 = self.v0[0]+t*self.dv[0]; xt1 = self.v0[1]+t*self.dv[1]
        mx = (self.v0[0]+self.v1[0])/2; my = (self.v0[1]+self.v1[1])/2
        Td = (self.dv[0]/self.L, self.dv[1]/self.L)
        ell = Td[0]*(xt0-mx)+Td[1]*(xt1-my)
        return (self.H+z)*mp.cos(u) + ell*mp.sin(u)
    def tent(self, u):
        return self.H*mp.cos(u) + (self.L/2)*abs(mp.sin(u))
    def Dval(self, s, tau):
        u = tau*s
        return mp.log(self.h(u, 1/(1-tau))/self.tent(u))/tau

def extract(mm, s_list, tau0):
    out = []
    def nev(pts):
        n = len(pts); Q = [[mp.mpf(0)]*n for _ in range(n)]
        for i in range(n): Q[i][0] = pts[i][1]
        for jj in range(1,n):
            for i in range(n-jj):
                xi, xj = pts[i][0], pts[i+jj][0]
                Q[i][jj] = (xj*Q[i][jj-1] - xi*Q[i+1][jj-1])/(xj-xi)
        return Q[0][n-1]
    Zo = Loc(mm)
    for s in s_list:
        vals = [Zo.Dval(s, tau0/2**j) for j in range(4)]
        taus = [tau0/2**j for j in range(4)]
        D1 = nev(list(zip(taus, vals)))
        D2 = nev([(t,(v-D1)/t) for t,v in zip(taus,vals)])
        D3 = nev([(t,(v-D1-D2*t)/t**2) for t,v in zip(taus,vals)])
        out.append((s,D1,D2,D3))
    return out

if __name__ == "__main__":
    mm = int(sys.argv[1])
    s_list = [mp.mpf(k)/4 for k in range(-96, 97)]
    prof = extract(mm, s_list, mp.mpf('0.006'))
    json.dump([[str(s),str(a),str(b),str(c)] for s,a,b,c in prof], open(f"ckpt_locprof_m{mm}.json","w"))
    for s,a,b,c in prof[::24]:
        print(f"m={mm}loc s={mp.nstr(s,6)}: D1={mp.nstr(a,16)} D2={mp.nstr(b,14)} D3={mp.nstr(c,14)}")
