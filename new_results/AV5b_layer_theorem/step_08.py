import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
with open("progress.log","a") as f:
    f.write("step_00 setup+reading: AV5 result.md, paper diamond-end (P'''(1)=6sqrt2 z3-2sqrt2 pi^2, pure-weight P/sqrt2=4+z2 t^2+z3 t^3-33/16 z4 t^4+...), C_m law read\n")

step = r'''
import mpmath as mp
mp.mp.dps = 50

def S_p(phi, m, p):
    return mp.nsum(lambda i: abs(mp.cos(phi - i*mp.pi/m))**p, [1, m]) if False else \
           mp.fsum(abs(mp.cos(phi - i*mp.pi/m))**p for i in range(1, m+1))

def rho_p(phi, m, p):
    return S_p(phi, m, p)**(-mp.mpf(1)/p)

def h_p(theta, m, p):
    # support = max_phi rho_p(phi) cos(theta - phi)
    f = lambda phi: rho_p(phi, m, p)*mp.cos(theta - phi)
    # reduce theta into the fundamental sector; use findroot on derivative
    def df(phi):
        return mp.diff(f, phi)
    best = None
    # coarse grid to bracket
    N = 720
    phis = [2*mp.pi*k/N for k in range(N)]
    vals = [f(ph) for ph in phis]
    k0 = max(range(N), key=lambda k: vals[k])
    a, b = phis[(k0-1) % N], phis[(k0+1) % N]
    try:
        root = mp.findroot(df, phis[k0], tol=mp.mpf('1e-40'))
        return f(root)
    except Exception:
        return vals[k0]

# polygon (p=1) support: tent data
for m in [2,3]:
    p = 1
    th0 = mp.pi/2 - mp.pi/(2*m)   # edge normal
    H = h_p(th0, m, 1)
    # slopes: derivative from each side
    du = mp.mpf('1e-20')
    sl_r = (h_p(th0+du, m, 1) - H)/du
    sl_l = (h_p(th0-du, m, 1) - H)/(-du)
    print(f"m={m}: th0={mp.nstr(th0,8)} H={mp.nstr(H,20)} sin(pi/2m)={mp.nstr(mp.sin(mp.pi/(2*m)),20)}")
    print(f"   slopes +side={mp.nstr(sl_r,10)} -side={mp.nstr(sl_l,10)}  T=H tan(pi/2m)={mp.nstr(H*mp.tan(mp.pi/(2*m)),10)}")

# overshoot at corner vs d
for m in [2,3,4]:
    th0 = mp.pi/2 - mp.pi/(2*m)
    H = h_p(th0, m, 1)
    for d in [mp.mpf('0.2'), mp.mpf('0.1'), mp.mpf('0.05')]:
        p = 1+d
        tau = 1 - 1/p
        r = h_p(th0, m, p)/H
        print(f"m={m} d={mp.nstr(d,4)} tau={mp.nstr(tau,6)}: h_p/H at corner = {mp.nstr(r,12)}  log(r)/tau = {mp.nstr(mp.log(r)/tau,12)}  (ln2={mp.nstr(mp.log(2),8)})")
'''
with open("step_02.py","w") as f: f.write(step)
exec(step)
