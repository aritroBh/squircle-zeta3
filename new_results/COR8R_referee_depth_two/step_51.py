import mpmath as mp
mp.mp.dps = 25
c0 = mp.mpf('0.31')
B = 2/(1+c0)
G2 = mp.log(1-mp.mpf(1)/2)
GB = mp.log(1-1/B)
G2B = mp.quad(lambda t1: mp.log(1-t1/B)/(t1-2), [mp.mpf(0),mp.mpf(1)])
GB2 = mp.quad(lambda t1: mp.log(1-t1/2)/(t1-B), [mp.mpf(0),mp.mpf(1)])
lhs = G2*GB
rhs = G2B+GB2
print("G(2)G(B) =", mp.nstr(lhs,20), " G(2,B)+G(B,2) =", mp.nstr(rhs,20), " diff:", mp.nstr(abs(lhs-rhs),5))
# Now check S(G(2,B)) via its DE in c: d/dlog c G(2,B;1) =? G(B)*dlog(letter1) + G(2)*dlog(letter2)
# letter1 = (1-2)/(2-B) = -(1+c)/(2c): dlog = dlog(1+c)-dlog c -> c/(1+c) - 1 = -1/(1+c)
# letter2 = (2-B)/B = c: dlog c -> 1
h = mp.mpf('1e-12')
def G2B_of(cv):
    Bv = 2/(1+cv)
    return mp.quad(lambda t1: mp.log(1-t1/Bv)/(t1-2), [mp.mpf(0),mp.mpf(1)])
dG = (G2B_of(c0*mp.exp(h))-G2B_of(c0*mp.exp(h).neg() if False else c0*mp.exp(-h)))/(2*h)
pred = GB*(-1/(1+c0)) + G2*(1)
print("d G(2,B)/dlog c:", mp.nstr(dG,20), " pred:", mp.nstr(pred,20), " diff:", mp.nstr(abs(dG-pred),5))
