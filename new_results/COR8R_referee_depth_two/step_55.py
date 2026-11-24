import mpmath as mp
mp.mp.dps = 25
def G3(w, cval):
    B = 2/(1+cval)
    w = [mp.mpf(2) if x==2 else (mp.mpf(0) if x==0 else B) for x in w]
    return mp.quad(lambda t1: mp.quad(lambda t2: mp.quad(lambda t3: 1/(t3-w[2]), [mp.mpf(0),t2])/(t2-w[1]), [mp.mpf(0),t1])/(t1-w[0]), [mp.mpf(0),mp.mpf(1)])
def G2f(w, cval):
    B = 2/(1+cval)
    w = [mp.mpf(2) if x==2 else (mp.mpf(0) if x==0 else B) for x in w]
    return mp.quad(lambda t1: mp.quad(lambda t2: 1/(t2-w[1]), [mp.mpf(0),t1])/(t1-w[0]), [mp.mpf(0),mp.mpf(1)])
c0 = mp.mpf('0.31'); h = mp.mpf('1e-11')
def dnum(f, w):
    return (f(w, c0*mp.exp(h))-f(w, c0*mp.exp(-h)))/(2*h)
B0 = 2/(1+c0)
G0B  = -mp.polylog(2, 1/B0)      # G(0,B;1) = -Li2(1/B)
G02  = -mp.polylog(2, mp.mpf('0.5'))  # G(0,2;1) = -Li2(1/2)
G2B  = G2f([2,B0], c0); GB2 = G2f([B0,2], c0)
print("check s20B: dG(2,0,B)/dlogc vs G(2,B)*c/(1+c)")
print("  ", mp.nstr(dnum(G3,[2,0,B0]),15), mp.nstr(G2B*c0/(1+c0),15))
print("check s02B: dG(0,2,B)/dlogc vs -G(0,B)/(1+c) + G(0,2)")
print("  ", mp.nstr(dnum(G3,[0,2,B0]),15), mp.nstr(-G0B/(1+c0)+G02,15))
print("check s0B2: dG(0,B,2)/dlogc vs G(B,2)*c/(1+c) - G(0,2) + G(0,B)/(1+c)")
print("  ", mp.nstr(dnum(G3,[0,B0,2]),15), mp.nstr(GB2*c0/(1+c0)-G02+G0B/(1+c0),15))
