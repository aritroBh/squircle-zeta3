import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
exec(open("step_27.py").read().split("for mm, full in")[0])
mp.mp.dps = 60
v0, dv, ns, idx, sig, al, be, b, B0, L, H = setup(3, False)
print("v0:", [mp.nstr(x,8) for x in v0], "dv:", [mp.nstr(x,8) for x in dv])
print("ns:", [[mp.nstr(x,6) for x in n] for n in ns])
print("idx:", idx, "sig:", sig)
print("al:", [mp.nstr(x,8) for x in al], "be:", [mp.nstr(x,8) for x in be])
print("b:", [mp.nstr(x,8) for x in b], "B0:", mp.nstr(B0,10), "L:", mp.nstr(L,10), "H:", mp.nstr(H,10))
z1p, z2p, L = z1z2(3, False)
for t in [mp.mpf('0.05'), mp.mpf('0.3'), mp.mpf('0.5')]:
    print(f"t={mp.nstr(t,4)}: z1p={mp.nstr(z1p(t),12)} z2p={mp.nstr(z2p(t),12)} prod={mp.nstr(z1p(t)*z2p(t),12)}")
I = mp.quad(lambda t: z1p(t)*z2p(t), [0, mp.mpf('0.5'), 1])
print("I =", mp.nstr(I, 20))
