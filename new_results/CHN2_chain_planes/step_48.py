
import mpmath as mp
mp.mp.dps = 70
def single_V(u2,u3,u4,u5, u, M=3000):
    # V(u) = sum_m (-u5 u)_m/m! u4 u/(u4 u + m) K_m(u) with tail Richardson
    S3 = u2+u3+u4
    def Km(m):
        return mp.gamma(1+u3*u)*mp.gamma(1+(u2+u4)*u+m)/mp.gamma(1+S3*u+m)
    def term(m):
        return mp.rf(-u5*u,m)/mp.factorial(m) * (u4*u)/(u4*u+m) * Km(m)
    Ms = [M, M*3//4, M//2, M*5//12, M//3]
    vals=[]
    for Mi in Ms:
        vals.append(mp.fsum(term(m) for m in range(Mi)))
    # tail exponent: term ~ (u5 u)*... m^(-u5 u -1) * u4u/m * m^(-u3 u) = C m^(-2 - (u3+u5)u)
    s_exp = 2 + (u3+u5)*u
    n=len(Ms)
    A = mp.matrix(n,n)
    for i,Mi in enumerate(Ms):
        A[i,0]=1
        for j in range(1,n):
            A[i,j]=Mi**(1-s_exp)*Mi**(-(j-1))
    return mp.lu_solve(A, mp.matrix(vals))[0]

def extract_c(u2,u3,u4,u5, Wmax=6):
    hs=[mp.mpf(1)/300, mp.mpf(1)/600, mp.mpf(1)/1200, mp.mpf(1)/2400, mp.mpf(1)/4800, mp.mpf(1)/9600][:Wmax-1]
    Vs=[single_V(u2,u3,u4,u5,h) for h in hs]
    n=len(hs)
    A=mp.matrix(n,n)
    for i,h in enumerate(hs):
        for j in range(n):
            A[i,j]=h**(j+2)
    sol=mp.lu_solve(A, mp.matrix([V-1 for V in Vs]))
    return sol

for pt in [(1,1,1,1),(2,1,1,1)]:
    sol = extract_c(*pt)
    print(pt, [mp.nstr(s,12) for s in sol])
    # compare with zeta combos
    print("  c3/z3:", mp.nstr(sol[1]/mp.zeta(3),12), " c4/z4:", mp.nstr(sol[2]/mp.zeta(4),12))
    print("  c5/z5:", mp.nstr(sol[3]/mp.zeta(5),12), " c5/(z2z3):", mp.nstr(sol[3]/(mp.zeta(2)*mp.zeta(3)),12))
