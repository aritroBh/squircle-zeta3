import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
exec(open("step_27.py").read().split("for mm, full in")[0])
mp.mp.dps = 60
v0, dv, ns, idx, sig, al, be, b, B0, L, H = setup(3, False)
def parts(t):
    out = []
    for i, s in zip(idx, sig):
        a = al[i] + be[i]*t
        aa = abs(a)
        Li = mp.log(aa) if aa > 0 else mp.mpf('-1e300')
        out.append((a, aa, Li, s, be[i], b[i]))
    return out
def z1(t):
    P = parts(t); return -mp.fsum(aa*Li for a, aa, Li, s, bi, bb in P)/B0
def z1p(t):
    P = parts(t); return -mp.fsum(bi*s*(Li+1) for a, aa, Li, s, bi, bb in P)/B0
t = mp.mpf('1e-6')
P = parts(t)
z = z1(t); zp = z1p(t)
print("z1:", mp.nstr(z,15), "z1p:", mp.nstr(zp,15))
for (a, aa, Li, s, bi, bb) in P:
    print(f"a={mp.nstr(a,6)} sig={s} be={mp.nstr(bi,4)} b={mp.nstr(bb,4)} Li={mp.nstr(Li,8)}")
    print(f"   T1a = sig*zp*b*(1+Li) = {mp.nstr(s*zp*bb*(1+Li),12)}")
    print(f"   T1b = sig*z*b*be/a = {mp.nstr(s*z*bb*bi/a,12)}")
    print(f"   T2  = sig*be*(Li^2/2+Li) = {mp.nstr(s*bi*(Li**2/2+Li),12)}")
