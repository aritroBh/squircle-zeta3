# step_73: direct vs inversion at the sols arguments
import mpmath as mp
mp.mp.dps = 80
def Li3(x): return mp.re(mp.polylog(3, x))
def Li2(x): return mp.re(mp.polylog(2, x))
def L3_direct(x):
    x = mp.mpf(x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(abs(1-x))
def L3_inv(x):
    x = mp.mpf(x)
    if x > 1: return L3_direct(1/x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(1-x)
for xx in ['2','3','4','5','6','7','8','9','10','15','16','21','25','28','36','10/9','5/4']:
    x = mp.mpf(xx)
    d, i = L3_direct(x), L3_inv(x)
    print(f"x={xx}: direct={mp.nstr(d,25)} inv={mp.nstr(i,25)} diff={mp.nstr(abs(d-i),4)}")
