import mpmath as mp
mp.mp.dps = 25
def G3(w, cval):
    B = 2/(1+cval)
    w = [mp.mpf(2) if x==2 else (mp.mpf(0) if x==0 else B) for x in w]
    return mp.quad(lambda t1: mp.quad(lambda t2: mp.quad(lambda t3: 1/(t3-w[2]), [mp.mpf(0),t2])/(t2-w[1]), [mp.mpf(0),t1])/(t1-w[0]), [mp.mpf(0),mp.mpf(1)])
for c0 in [mp.mpf('0.31'), mp.mpf('0.47'), mp.mpf('0.62')]:
    B0 = 2/(1+c0)
    G2 = mp.log(mp.mpf('0.5'))
    GB0 = mp.polylog(2,1/B0)
    lhs = G2*GB0
    rhs = -(G3([2,0,B0],c0)+G3([0,2,B0],c0)+G3([0,B0,2],c0))
    print(f"c={c0}: lhs={mp.nstr(lhs,15)} rhs={mp.nstr(rhs,15)} diff={mp.nstr(abs(lhs-rhs),3)}")
