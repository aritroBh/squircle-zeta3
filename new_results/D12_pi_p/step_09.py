
import mpmath as mp
import pickle

mp.mp.dps = 110
D = pickle.load(open('/tmp/Nvals.pkl','rb'))
INT = D['INT']
log2 = mp.log(2)
KMAX = 44
# true N_k = [h^(k+1)] (pi_p/4) = (-log2)^(k+1)/(k+1)! + sum_{m=0}^{k-1} INT[m]
N = [mp.mpf(0)]*(KMAX+2)
for k in range(0,KMAX+1):
    N[k] = (-log2)**(k+1)/mp.factorial(k+1) + sum(INT[0:k])  # INT[0..k-1]

# cross-check first 10 coefficients by direct fit of pi_{1/h}
mp.mp.dps = 150
def pi_cosh_h(h):
    f = lambda u: mp.cosh(u/2)**(-1-h) * mp.cosh(u*(1-h)/(2*h))**h
    return 2*h*mp.quad(f, [0, mp.inf])
# evaluate at h = 1/p for p = P0..P0+11, fit
import numpy as np
ps = [60+3*i for i in range(14)]
fs = [pi_cosh_h(1/mp.mpf(p))/4 for p in ps]
# finite-difference style: subtract known terms sequentially using mpmath matrix solve
M = len(ps)
A = mp.matrix(M,11)
for a,p in enumerate(ps):
    h = 1/mp.mpf(p)
    for k in range(11):
        A[a,k] = h**(k+1)
b = mp.matrix(M,1)
for a in range(M): b[a] = fs[a]-1
coef = mp.lu_solve(A, b)
mp.mp.dps = 40
print("k : series-arithmetic N_k        vs  direct-fit")
for k in range(11):
    print(f"N_{k}: {mp.nstr(N[k],25)}   {mp.nstr(coef[k],25)}   diff={mp.nstr(abs(N[k]-coef[k]),2)}")
pickle.dump({'N':N}, open('/tmp/Nfinal.pkl','wb'))
