
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
    L = mp.sqrt((v1[0]-v0[0])**2 + (v1[1]-v0[1])**2)
    ns = [(mp.cos(i*mp.pi/m), mp.sin(i*mp.pi/m)) for i in range(1, m+1)]
    return rho, v0, v1, th0, N, H, L, ns

def make_z(m, full):
    rho, v0, v1, th0, N, H, L, ns = geom(m)
    j0, j1 = m-1, m-2
    def ab(t):
        xt = (v0[0]+t*(v1[0]-v0[0]), v0[1]+t*(v1[1]-v0[1]))
        a = [n[0]*xt[0]+n[1]*xt[1] for n in ns]
        b = [n[0]*N[0]+n[1]*N[1] for n in ns]
        return a, b
    def F(t, z, p):
        a, b = ab(t)
        s = abs(a[j0]+z*b[j0])**p + abs(a[j1]+z*b[j1])**p
        for i in range(m):
            if i in (j0, j1): continue
            s += abs(a[i]+z*b[i])**p if full else abs(a[i])
        return s - 1
    def Fz(t, z, p):
        a, b = ab(t)
        s = b[j0]*mp.sign(a[j0]+z*b[j0])*abs(a[j0]+z*b[j0])**(p-1) if a[j0]+z*b[j0] != 0 else 0
        s += b[j1]*mp.sign(a[j1]+z*b[j1])*abs(a[j1]+z*b[j1])**(p-1) if a[j1]+z*b[j1] != 0 else 0
        s *= p
        if full:
            for i in range(m):
                if i in (j0,j1): continue
                aa = a[i]+z*b[i]
                if aa != 0:
                    s += p*b[i]*mp.sign(aa)*abs(aa)**(p-1)
        return s
    def z_solve(t, p):
        if t <= 0 or t >= 1: return mp.mpf(0)
        # Newton from z1 estimate
        a, b = ab(t)
        S1 = sum(abs(ai)*mp.log(abs(ai)) if ai != 0 else mp.mpf(0) for ai in a)
        Sb = sum(mp.sign(ai)*bi for ai, bi in zip(a, b))
        if not full:
            # frozen-rest derivative balance: Sb only over the two local kinks
            Sb = mp.sign(a[j0])*b[j0] + mp.sign(a[j1])*b[j1]
        z = -mp.mpf(p-1)*S1/Sb if Sb != 0 else mp.mpf(0)
        for _ in range(60):
            fz = F(t, z, p)
            if abs(fz) < mp.mpf('1e-45'): break
            dz = fz/Fz(t, z, p)
            z -= dz
            if abs(dz) < mp.mpf('1e-45'): break
        return z
    def zp(t, p):
        # z' = -F_t/F_z with F_t via complex-free analytic diff of a_i in t: a_i' = <n_i, v1-v0>
        z = z_solve(t, p)
        a, b = ab(t)
        ap = [n[0]*(v1[0]-v0[0]) + n[1]*(v1[1]-v0[1]) for n in ns]
        Ft = mp.mpf(0)
        for i in range(m):
            if not full and i not in (j0, j1):
                Ft += mp.sign(a[i])*ap[i]
                continue
            aa = a[i] + z*b[i]
            if aa != 0:
                Ft += p*mp.sign(aa)*abs(aa)**(p-1)*ap[i]
        return z, -Ft/Fz(t, z, p)
    return z_solve, zp, L

def dP_total(m, d, full):
    z_solve, zp, L = make_z(m, full)
    p = 1+d
    f = lambda t: mp.sqrt(L**2 + zp(t, p)[1]**2) - L
    return 2*m*mp.quad(f, [0, mp.mpf('0.5'), 1])

if __name__ == "__main__":
    mm = int(sys.argv[1]); full = int(sys.argv[2])
    ds = [mp.mpf(x) for x in sys.argv[3].split(",")]
    res = []
    for d in ds:
        v = dP_total(mm, d, full)
        res.append((d, v))
        print(f"m={mm} full={full} d={mp.nstr(d,4)}: DeltaP={mp.nstr(v,35)}", flush=True)
    json.dump([[str(d), str(v)] for d, v in res], open(f"ckpt_arc_m{mm}_f{full}.json","w"))
