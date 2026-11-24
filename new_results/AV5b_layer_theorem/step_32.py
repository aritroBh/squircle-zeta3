import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
exec(open("step_27.py").read().split("for mm, full in")[0])

z1p, z2p, L = z1z2(2, False)
for t in [mp.mpf('0.1'), mp.mpf('0.5'), mp.mpf('0.9'), mp.mpf('1e-6')]:
    print(f"t={mp.nstr(t,6)}: z1p={mp.nstr(z1p(t),12)} z2p={mp.nstr(z2p(t),12)}")
# finite-difference check of z2p
z1, z2 = None, None
v0, dv, ns, idx, sig, al, be, b, B0, L, H = setup(2, False)
print("idx:", idx, "sig:", sig, "B0:", mp.nstr(B0,10))
def z1f(t):
    P = [(al[i]+be[i]*t, abs(al[i]+be[i]*t), 0, sig[k], be[i], b[i]) for k,i in enumerate(idx)]
    N1 = mp.fsum(aa*mp.log(aa) for a, aa, *_ in P)
    return -N1/B0
print("z1(0.5) fd check:", mp.nstr(z1f(mp.mpf('0.5')),12))
print("z1p(0.5):", mp.nstr(z1p(mp.mpf('0.5')),12), " fd:", mp.nstr(mp.diff(z1f, mp.mpf('0.5')),12))
