
import pickle
import mpmath as mp
mp.mp.dps = 250
acc = pickle.load(open("local_block_exact.pkl","rb"))
r5 = mp.sqrt(5); phi = (1+r5)/2
l2 = mp.log(2); lf = mp.log(phi)
z2 = mp.pi**2/6; z3 = mp.zeta(3); z5 = mp.zeta(5)
ZVm = {'1': mp.mpf(1), 'z2': z2, 'z3': z3, 'z2^2': z2**2, 'z5': z5, 'z2z3': z2*z3}
def q5mpf(z):
    return mp.mpf(z.p.numerator)/mp.mpf(z.p.denominator) + (mp.mpf(z.q.numerator)/mp.mpf(z.q.denominator))*r5
val = sum(q5mpf(acc[k])*ZVm[k[0]]*l2**k[1]*lf**k[2] for k in acc)
local250 = mp.mpf(pickle.load(open("local250.pkl","rb")))
print("exact  :", mp.nstr(val, 55))
print("quad   :", mp.nstr(local250, 55))
print("residual:", mp.nstr(abs(val-local250), 3))
# cross-check with 110-digit itab block value
prev = mp.mpf(pickle.load(open("w5L_num.pkl","rb"))['blocks']['local'])
print("vs 110-digit itab block:", mp.nstr(abs(val-prev), 3))
