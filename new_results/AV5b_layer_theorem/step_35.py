import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
exec(open("step_27.py").read().split("for mm, full in")[0])
mp.mp.dps = 60
z1p, z2p, L = z1z2(3, False)
for t in [mp.mpf('1e-3'), mp.mpf('1e-6'), mp.mpf('1e-10'), mp.mpf('1e-20')]:
    print(f"t={mp.nstr(t,4)}: z1p={mp.nstr(z1p(t),12)} z2p={mp.nstr(z2p(t),12)} prod={mp.nstr(z1p(t)*z2p(t),12)}")
I1 = mp.quad(lambda t: z1p(t)*z2p(t), [mp.mpf('0.001'), mp.mpf('0.5'), mp.mpf('0.999')])
I2 = mp.quad(lambda t: z1p(t)*z2p(t), [mp.mpf('1e-8'), mp.mpf('0.001')])
print("I(0.001..0.999) =", mp.nstr(I1, 15))
print("I(1e-8..0.001) =", mp.nstr(I2, 15))
