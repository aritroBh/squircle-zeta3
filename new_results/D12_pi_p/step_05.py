# step_05.py -- TASK 3: assemble N_k (k<=43) from PSI via J-integrals at 110 digits; cross-check by direct fit
import mpmath as mp, pickle
mp.mp.dps = 110
PSI = pickle.load(open('/tmp/psi.pkl','rb'))['PSI']; KMAX=44
pairs=set()
for co in PSI:
    for (i,j,k) in co: pairs.add((j,i))
J={}
for (j,i) in sorted(pairs):
    J[(j,i)] = mp.quad(lambda s: s**i*mp.log(1+mp.e**(-s))**j, [0, mp.inf])
log2=mp.log(2)
INT=[mp.mpf(0)]*(KMAX+1)
for m in range(KMAX+1):
    INT[m] = sum(mp.mpf(v.numerator)/mp.mpf(v.denominator)*J[(j,i)]*log2**k
                 for (i,j,k),v in PSI[m].items())
N=[mp.mpf(0)]*(KMAX+2)
for k in range(0,KMAX+1):
    N[k] = (-log2)**(k+1)/mp.factorial(k+1) + sum(INT[0:k])
for k in range(12): print(f"N_{k:2d} = {mp.nstr(N[k],40)}")
pickle.dump({'N':N,'INT':INT,'J':J}, open('/tmp/Nvals.pkl','wb'))
# cross-check first 11 coefficients by direct fit of pi_{1/h}/4 at 150 digits
mp.mp.dps=150
def pi_h(h):
    f = lambda u: mp.cosh(u/2)**(-1-h)*mp.cosh(u*(1-h)/(2*h))**h
    return 2*h*mp.quad(f,[0,mp.inf])/4
ps=[60+3*i for i in range(14)]; fs=[pi_h(1/mp.mpf(p)) for p in ps]
A=mp.matrix(14,11); b=mp.matrix(14,1)
for a,p in enumerate(ps):
    h=1/mp.mpf(p)
    for k in range(11): A[a,k]=h**(k+1)
    b[a]=fs[a]-1
coef=mp.lu_solve(A,b)
mp.mp.dps=40
print("direct fit vs series arithmetic:", [mp.nstr(abs(coef[k]-N[k]),2) for k in range(11)])
