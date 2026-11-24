
import mpmath as mp
import pickle
mp.mp.dps = 110
N = pickle.load(open('/tmp/Nfinal.pkl','rb'))['N']
pi=mp.pi

def pi_cosh(p):
    p=mp.mpf(p); h=1/p
    f = lambda u: mp.cosh(u/2)**(-1-h)*mp.cosh(u*(1-h)/(2*h))**h
    return 2*h*mp.quad(f,[0,mp.inf])

print("=== Optimal truncation of the asymptotic series at p=10 (h=0.1) ===")
true = pi_cosh(10)/4
h = mp.mpf('0.1')
s = mp.mpf(1); best=(0,mp.mpf(1))
for k in range(0,42):
    s += N[k]*h**(k+1)
    err = abs(s-true)
    if err < best[1]: best=(k,err)
    if k in [3,5,10,20,30,31,32,35,40]:
        print(f"  trunc k={k:2d}: err={mp.nstr(err,3)}")
print(f"  best k={best[0]}, err={mp.nstr(best[1],3)};  pi/h = {mp.nstr(pi/h,6)}, e^(-pi/h) = {mp.nstr(mp.e**(-pi/h),3)}")
print(f"  pi_10 = {mp.nstr(4*true,60)}")
