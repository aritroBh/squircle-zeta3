from mpmath import *
mp.dps = 90

def consts(th):
    z = sin(th/2)**2; c = cos(th); w = (1-z)/z; l = log(1-z); u = w*l
    S12 = zeta(3) - polylog(3,1-z) + l*polylog(2,1-z) + l**2*log(z)/2
    a1 = -1-u
    a2 = zeta(2)+mpf(1)/2-polylog(2,z)-u**2/2
    a3 = 2*w*S12 - (1+u**3)/3
    return z,c,a1,a2,a3

def R(th,e):
    z = sin(th/2)**2
    return pi*e/sin(pi*e)*2*z*hyp2f1(-e,1+e,2,z)

def Lexact(th,n):
    n = mpf(n); e = 1/(n-1); c = cos(th)
    s = mpf(0); m = 1; prev_small = False
    while True:
        Lam = (legendre(m-1,c)-legendre(m+1,c))/(2*m+1)
        t = (-1)**m*e/(m-e)*Lam*n**(-m)
        s += t
        small = abs(t) < mpf(10)**(-mp.dps+5)
        if small and prev_small and m>3: break   # two consecutive tiny terms (odd/even symmetry at pi/2)
        prev_small = small
        m += 1
    return 1+c+R(th,e)*n**(-e)+s

def Ldirect(th,n):
    c = cos(th)
    return quad(lambda x: sqrt(1+2*c*n*x**(n-1)+(n*x**(n-1))**2),linspace(0,1,12))

def a_cauchy(th,N=128,r=mpf(1)/4):
    # Taylor coefficients of log(R/2z) via Cauchy trapezoid
    z = sin(th/2)**2
    def lg(e): return log(pi*e/sin(pi*e)*hyp2f1(-e,1+e,2,z))
    vals=[lg(r*exp(2j*pi*k/N)) for k in range(N)]
    return [ (sum(vals[k]*exp(-2j*pi*k*j/N) for k in range(N))/N/r**j).real for j in range(1,4)]

if __name__=="__main__":
    print("== formula vs literal arc-length quadrature ==")
    for th in (pi/3,pi/2,2*pi/3):
        for n in (3,10):
            print(nstr(th,5),n,nstr(Lexact(th,n)-Ldirect(th,n),3))
    print("== Taylor coefficients: closed form vs Cauchy ==")
    for th in (pi/3,pi/2,2*pi/3):
        z,c,a1,a2,a3 = consts(th); ac = a_cauchy(th)
        print(nstr(th,5),[nstr(abs(x-y),3) for x,y in zip((a1,a2,a3),ac)], "a1,a2,a3 =",nstr(a1,10),nstr(a2,10),nstr(a3,10))
    print("== table: D = (L - [1+c+2z(mu/n)^eps + c2/n^2]) n^3/log n ;  pred = -2 z a2 + d/log n ;  r4 = (E n^3 - (d-2 z a2 log n)) n/log^2 n ==")
    for th,name in ((pi/3,"pi/3"),(pi/2,"pi/2"),(2*pi/3,"2pi/3"),(mpf('0.1'),"0.1"),(pi-mpf('0.1'),"pi-0.1")):
        z,c,a1,a2,a3 = consts(th)
        mu = exp(a1); c2 = 2*z*(a2-1+z)
        d = 2*z*(a2*(2+a1)+a3) - z*(3-z-2*z**2)
        print("theta =",name," z=",nstr(z,8)," -2za2=",nstr(-2*z*a2,10))
        print("%8s %16s %16s %14s %12s"%("n","D","pred","D/z","r4"))
        for k in range(2,7):
            n = mpf(10)**k; e = 1/(n-1)
            E = Lexact(th,n) - (1+c+2*z*(mu/n)**e + c2/n**2)
            ln = log(n)
            D = E*n**3/ln
            pred = -2*z*a2 + d/ln
            r4 = (E*n**3-(d-2*z*a2*ln))*n/ln**2
            print("%8s %16s %16s %14s %12s"%("1e%d"%k,nstr(D,10),nstr(pred,10),nstr(D/z,10),nstr(r4,6)))
