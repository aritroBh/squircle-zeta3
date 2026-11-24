import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
step = r'''
import mpmath as mp
mp.mp.dps = 60

def setup(m, full):
    rho = mp.tan(mp.pi/(2*m))
    v0 = (mp.mpf(0), rho)
    v1 = (rho*mp.sin(mp.pi/m), rho*mp.cos(mp.pi/m))
    th0 = mp.pi/2 - mp.pi/(2*m)
    N = (mp.cos(th0), mp.sin(th0))
    H = mp.sin(mp.pi/(2*m))
    dv = (v1[0]-v0[0], v1[1]-v0[1])
    L = mp.sqrt(dv[0]**2+dv[1]**2)
    ns = [(mp.cos(i*mp.pi/m), mp.sin(i*mp.pi/m)) for i in range(1, m+1)]
    jset = {m-1, m-2} if not full else set(range(m))
    al = [n[0]*v0[0]+n[1]*v0[1] for n in ns]        # a_i(0)
    be = [n[0]*dv[0]+n[1]*dv[1] for n in ns]        # a_i'
    b  = [n[0]*N[0]+n[1]*N[1] for n in ns]
    idx = sorted(jset)
    sig = [mp.sign(al[i]) if al[i] != 0 else mp.sign(al[i]+be[i]*mp.mpf('0.3')) for i in idx]
    B0 = mp.fsum(s*bb for s, bb in zip(sig, [b[i] for i in idx]))
    return v0, dv, ns, idx, sig, al, be, b, B0, L, H

def z1z2(m, full):
    v0, dv, ns, idx, sig, al, be, b, B0, L, H = setup(m, full)
    def parts(t):
        out = []
        for i, s in zip(idx, sig):
            a = al[i] + be[i]*t
            aa = abs(a)
            Li = mp.log(aa) if aa > 0 else mp.mpf('-1e300')
            out.append((a, aa, Li, s, be[i], b[i]))
        return out
    def z1(t):
        P = parts(t)
        N1 = mp.fsum(aa*Li for a, aa, Li, s, bi, bb in P)
        return -N1/B0
    def z1p(t):
        P = parts(t)
        return -mp.fsum(bi*s*(Li+1) for a, aa, Li, s, bi, bb in P)/B0
    def z2(t):
        P = parts(t)
        z = z1(t)
        N2 = mp.fsum(s*(z*bb*(1+Li)) + aa*Li**2/2 for a, aa, Li, s, bi, bb in P)
        return -N2/B0
    def z2p(t):
        P = parts(t)
        z = z1(t); zp = z1p(t)
        T1 = mp.fsum(s*(zp*bb*(1+Li) + z*bb*bi/a) for a, aa, Li, s, bi, bb in P)
        T2 = mp.fsum(s*bi*(Li**2/2 + Li) for a, aa, Li, s, bi, bb in P)
        return -(T1+T2)/B0
    return z1p, z2p, L

def A3_model(m, full):
    z1p, z2p, L = z1z2(m, full)
    f = lambda t: z1p(t)*z2p(t)
    I = mp.quad(f, [0, mp.mpf('0.25'), mp.mpf('0.5'), mp.mpf('0.75'), 1])
    return 2*m/L*I, I

for mm, full in [(2,0),(3,0),(4,0)]:
    A3, I = A3_model(mm, full)
    print(f"m={mm} full={full}: A3 = {mp.nstr(A3, 30)}   6A3 = {mp.nstr(6*A3, 30)}")
    if mm == 2:
        print("   target sqrt2(z3-pi^2/3):", mp.nstr(mp.sqrt(2)*(mp.zeta(3)-mp.pi**2/3), 30))
'''
with open("step_27.py","w") as f: f.write(step)
exec(step)
