from cor8r_num import *
setdps(30)
c0 = mp.mpf('0.31'); B0 = 2/(1+c0)
def G3(w, cval):
    B = 2/(1+cval)
    w = [mp.mpf(2) if x==2 else (mp.mpf(0) if x==0 else B) for x in w]
    return mp.quad(lambda t1: mp.quad(lambda t2: mp.quad(lambda t3: 1/(t3-w[2]), [mp.mpf(0),t2])/(t2-w[1]), [mp.mpf(0),t1])/(t1-w[0]), [mp.mpf(0),mp.mpf(1)])
print("G(2,0,B;1) =", mp.nstr(G3([2,0,B0],c0),18), " Li_{2,1}(2/B,1/2) =", mp.nstr(Li_ab(2,1,2/B0,mp.mpf('0.5')),18))
print("G(0,2,B;1) =", mp.nstr(G3([0,2,B0],c0),18), " Li_{1,2}(2/B,1/2) =", mp.nstr(Li_ab(1,2,2/B0,mp.mpf('0.5')),18))
print("G(0,B,2;1) =", mp.nstr(G3([0,B0,2],c0),18), " Li_{1,2}(B/2,1/B) =", mp.nstr(Li_ab(1,2,B0/2,1/B0),18))
print("G(B,0,2;1) =", mp.nstr(G3([B0,0,2],c0),18), " Li_{2,1}(B/2,1/B) =", mp.nstr(Li_ab(2,1,B0/2,1/B0),18))
