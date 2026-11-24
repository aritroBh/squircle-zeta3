
import mpmath as mp, time, pickle
mp.mp.dps = 70
def volA_complex(u1,u2,u3,u4,u5, M=2500):
    def one(a1,a2,a4,a5):
        S4 = a1+u3+a5+3
        m0 = a2*mp.beta(a2+a4, 1+a1+u3)
        def alpha(m):
            return mp.rf(-a5,m)/mp.factorial(m)*(a4)/(a4+m)*a2*mp.beta(a2+a4+m, 1+a1+u3)
        Ms=[M, M*4//5, M*3//5, M//2, M*2//5]
        vals=[mp.fsum(alpha(m) for m in range(1,Mi)) for Mi in Ms]
        n=len(Ms)
        Am=mp.matrix(n,n)
        for i,Mi in enumerate(Ms):
            Am[i,0]=1
            for j in range(1,n): Am[i,j]=Mi**(1-S4)*Mi**(-(j-1))
        return m0 + mp.lu_solve(Am, mp.matrix(vals))[0]
    return one(u1,u2,u4,u5)+one(u5,u4,u2,u1)

a=(1,2,4,3,5)
r = mp.mpf(1)/30
N=64
vals=[]
t0=time.time()
for j in range(N):
    th = 2*mp.pi*j/N
    uj = r*mp.e**(1j*th)
    vals.append(volA_complex(*[mp.mpc(x)*uj for x in a]))
    if j%16==0: print(f"{j}/{N} {time.time()-t0:.0f}s", flush=True)
# FFT
coeffs=[]
for W in range(9):
    c = mp.fsum(vals[j]*mp.e**(-2*mp.pi*1j*W*j/N) for j in range(N))/N
    coeffs.append(c/ r**W)
pickle.dump([str(c) for c in coeffs], open('fft_check.pkl','wb'))
print("c0 =", mp.nstr(coeffs[0], 20))
for W in range(2,9):
    print(f"c{W} = {mp.nstr(coeffs[W], 25)}")
