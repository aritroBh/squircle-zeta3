
import mpmath as mp
mp.mp.dps = 50

def geom(m):
    rho = mp.tan(mp.pi/(2*m))
    v0 = (mp.mpf(0), rho)
    v1 = (rho*mp.sin(mp.pi/m), rho*mp.cos(mp.pi/m))
    th0 = mp.pi/2 - mp.pi/(2*m)
    N = (mp.cos(th0), mp.sin(th0))
    H = mp.sin(mp.pi/(2*m))
    L = mp.sqrt((v1[0]-v0[0])**2 + (v1[1]-v0[1])**2)
    Tdir = ((v1[0]-v0[0])/L, (v1[1]-v0[1])/L)
    ns = [(mp.cos(i*mp.pi/m), mp.sin(i*mp.pi/m)) for i in range(1, m+1)]
    return rho, v0, v1, th0, N, H, L, Tdir, ns

def x_t(t, v0, v1):
    return (v0[0]+t*(v1[0]-v0[0]), v0[1]+t*(v1[1]-v0[1]))

def z_solve(t, m, d, full, cache={}):
    # solve  |a_j0+zb_j0|^p + |a_j1+zb_j1|^p + rest = 1
    rho, v0, v1, th0, N, H, L, Tdir, ns = geom(m)
    p = 1+d
    xt = x_t(t, v0, v1)
    a = [n[0]*xt[0]+n[1]*xt[1] for n in ns]
    b = [n[0]*N[0]+n[1]*N[1] for n in ns]
    j0, j1 = m-1, m-2   # indices into ns: n_m -> ns[m-1], n_{m-1} -> ns[m-2]
    def constraint(z):
        s = abs(a[j0]+z*b[j0])**p + abs(a[j1]+z*b[j1])**p
        for i in range(m):
            if i in (j0, j1): continue
            s += abs(a[i])**p if full else abs(a[i])
        return s - 1
    c0 = constraint(mp.mpf(0))
    if abs(c0) < mp.mpf('1e-40'):
        return mp.mpf(0)
    if c0 < 0:
        z = mp.findroot(constraint, (mp.mpf(0), 3*H), solver='bisect', tol=mp.mpf('1e-42'))
    else:
        z = mp.findroot(constraint, (-2*H, mp.mpf(0)), solver='bisect', tol=mp.mpf('1e-42'))
    return z

def h_loc(theta, m, d, full):
    rho, v0, v1, th0, N, H, L, Tdir, ns = geom(m)
    u = theta - th0
    def val(t):
        z = z_solve(t, m, d, full)
        xt = x_t(t, v0, v1)
        ell = Tdir[0]*(xt[0]-(v0[0]+v1[0])/2) + Tdir[1]*(xt[1]-(v0[1]+v1[1])/2)
        return (H+z)*mp.cos(u) + ell*mp.sin(u)
    # coarse grid then refine
    Ngrid = 120
    ts = [mp.mpf(k)/Ngrid for k in range(1, Ngrid)]
    vals = [val(t) for t in ts]
    k0 = max(range(len(ts)), key=lambda k: vals[k])
    t0 = ts[k0]
    try:
        tr = mp.findroot(lambda t: mp.diff(val, t), t0, tol=mp.mpf('1e-40'))
        return val(tr)
    except Exception:
        return vals[k0]

def tent(u, m):
    rho, v0, v1, th0, N, H, L, Tdir, ns = geom(m)
    return H*mp.cos(u) + (L/2)*abs(mp.sin(u))

def dP_corner(m, d, full):
    W = mp.pi/(2*m)
    f = lambda u: h_loc(mp.pi/2 - mp.pi/(2*m) + u, m, d, full) - tent(u, m)
    return mp.quad(f, [-W, 0, W])

import sys
mm = int(sys.argv[1]) if len(sys.argv)>1 else 2
full = int(sys.argv[2]) if len(sys.argv)>2 else 0
ds = [mp.mpf('0.05'), mp.mpf('0.04'), mp.mpf('0.03'), mp.mpf('0.02'), mp.mpf('0.015')]
res = []
for d in ds:
    v = 2*mm*dP_corner(mm, d, full)
    res.append((d, v))
    print(f"m={mm} full={full} d={mp.nstr(d,4)}: DeltaP={mp.nstr(v,30)}", flush=True)
import json
json.dump([[str(d), str(v)] for d, v in res], open(f"ckpt_model_m{mm}_f{full}.json","w"))
