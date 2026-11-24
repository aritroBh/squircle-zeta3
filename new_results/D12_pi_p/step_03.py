# step_03.py -- TASK 2: closed forms N0..N3 (PROVED in writeup) and numerical confirmation
import mpmath as mp
mp.mp.dps = 80
pi = mp.pi; L2 = mp.log(2); z3 = mp.zeta(3)
N0 = -L2
N1 = L2**2/2
N2 = pi**2/24 - L2**3/6
N3 = L2**4/24 + pi**2/24*(1-L2) + z3/4

def pi_cosh(p):
    p = mp.mpf(p); h = 1/p
    f = lambda u: mp.cosh(u/2)**(-1-h)*mp.cosh(u*(1-h)/(2*h))**h
    return 2*h*mp.quad(f, [0, mp.inf])

print("N0..N3 =", [mp.nstr(x,40) for x in [N0,N1,N2,N3]])
for p in [40, 80, 160]:
    h = 1/mp.mpf(p); f = pi_cosh(p)/4
    S3 = 1 + N0*h + N1*h**2 + N2*h**3
    S4 = S3 + N3*h**4
    print(f"p={p}: f-S3={mp.nstr(f-S3,3)} (O(h^5), h^5={mp.nstr(h**5,3)})   f-S4={mp.nstr(f-S4,3)}")
