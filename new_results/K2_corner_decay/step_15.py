
import mpmath as mp
mp.mp.dps = 100
P=1024; R=mp.mpf('0.9')

# DFT self-test on exact data: f_p = u_p^7  => a_7 = R^7
def dft_test():
    half=P//2
    om=mp.exp(-2*mp.pi*mp.j/P)
    f=[(R*mp.exp(2*mp.pi*mp.j*p/P))**7 for p in range(half+1)]
    out=[]
    for k in [5,7,9]:
        s=mp.mpc(0); wk=om**k; wp=mp.mpc(1)
        for p in range(1,half):
            wp*=wk; s+=f[p]*wp
        a=(f[0]+((-1)**k)*f[half]+2*s.real)/P
        out.append((k, a*R**(-k), 2*s.imag/P))
    for k,a,im in out: print("k=%d  a=%s   ImErr=%.1e" % (k, mp.nstr(a,8), abs(im)))
dft_test()

# so DFT with wp*=wk accumulation: error check
# wp = om^(k p) via repeated mult: 511 muls each rel err 1e-100 -> fine.
# But wait: for k up to 300, wp accumulated 511 times: error ~ 511*1e-100 fine.
print("---- expected a_7 =", mp.nstr(R**7,8))
