from mpmath import mp, mpf, pi, sin, beta, gamma, log, matrix, lu_solve, factorial, nstr
mp.dps = 50
def pi_p(p):
    p = mpf(p); return 2*pi/(p*sin(pi/p))
def Nfun(p):
    p = mpf(p)
    return 4/(p*pi_p(p))*beta(3/p, 1-1/p)
hs = [mpf('0.01'), mpf('0.02')]
for h in hs:
    offs=list(range(-3,4))
    A = matrix(7,7); b = matrix(7,1)
    for i,j in enumerate(offs):
        for k in range(7): A[i,k]=(j*h)**k
        b[i]=Nfun(2+j*h)
    c = lu_solve(A,b)
    print("h=",h,"N''':", nstr(c[3]*6,35), " N'':", nstr(c[2]*2,35))
print("sympy N''': -5.35358499176667775189475449684")
print("sympy N'' : 0.457179055798648542005925120260")
